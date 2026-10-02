"""Re-implementation of the models in Kato (2022), PLOS ONE 17(9): e0275113.

"Islamic and capitalist economies: Comparison using econophysics models of
wealth exchange and redistribution."

Equation numbers below refer to the paper. Everything the paper leaves
unspecified is exposed as an explicit option (see `Params`) rather than being
decided silently; docs/kato2022_replication.md lists each one.

Exchange models (one random pair per time step):
    R  Eq (1)  Chatterjee et al. random exchange with savings rate lambda
    L  Eq (2)  loan interest: borrower i pays interest and bears all profit/loss
    J  Eq (3)  joint venture: both partners share the profit/loss rate delta
Redistribution (applied every `tp` steps, to all agents at once):
    T  Eq (4)  every agent gives xi*m_i, split equally among the other N-1
    Q  Eq (5-6) agents above m_max/q give xi*m_i to agents at or below m_max/q
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numba import njit

MODEL_R, MODEL_L, MODEL_J = 0, 1, 2
REDIST_NONE, REDIST_T, REDIST_Q = 0, 1, 2

_MODEL_CODES = {"R": MODEL_R, "L": MODEL_L, "J": MODEL_J}
_REDIST_CODES = {None: REDIST_NONE, "T": REDIST_T, "Q": REDIST_Q}
_INSOLV_CODES = {"none": 0, "floor": 1, "skip": 2, "cap": 3, "bankrupt": 4}


@dataclass(frozen=True)
class Params:
    """Parameters of one run. Defaults are the values used in Kato's figures."""

    model: str = "R"                 # "R", "L" or "J"
    n_agents: int = 1000             # N
    savings: float = 0.25            # lambda
    interest: float = 0.05           # rho (L model only)
    delta_width: float = 0.1         # delta_w: delta ~ U[-delta_w, +delta_w] (L, J)
    redistribution: str | None = None  # None, "T" or "Q"
    transfer_rate: float = 0.5       # xi
    transfer_period: int = 10**5     # tp
    quantile: float = 4.0            # q (Q model): threshold is m_max / q
    # Choices the paper does not pin down:
    shared_delta: bool = True        # J: one delta for the pair (Eq 3 uses one symbol)
    # Bankruptcy options (L model, insolvency="bankrupt" only; defaults = plain
    # liquidation in which the lender seizes everything the borrower has left).
    # At a default the borrower's assets are X = lam m_i + (1-lam)(1+d)(m_i+m_j) and
    # the lender's claim is C = (1-lam)(1+rho) m_j > X.
    #   bankruptcy_cost c: a fraction c of X is destroyed (legal fees, fire sale)
    #   recovery r:        the lender settles for min(r*C, (1-c)X); the borrower keeps
    #                      the rest. r<1 is a haircut / debt relief (Q 2:280), r=1 plain
    #                      liquidation. Wealth is conserved unless c>0.
    #   exclusion_period E: a defaulter cannot borrow for the next E steps (it can lend)
    recovery: float = 1.0
    bankruptcy_cost: float = 0.0
    exclusion_period: int = 0
    # L only: what happens when Eq (2) would push the borrower below zero.
    #   "bankrupt" borrower defaults: reset to 0, lender absorbs the shortfall
    #           (default; total wealth evolves exactly as under Eq 2; Result.defaults counts events)
    #   "skip"  trade rejected, both keep their wealth (inflates total wealth)
    #   "floor" borrower clamped at 0     "cap"  borrower pays only what it has
    #   "none"  Eq (2) literally: diverges to wealth ~ -1e18 once rho > 0
    insolvency: str = "bankrupt"
    samples: tuple[int, ...] = field(default_factory=tuple)  # times to record

    def __post_init__(self):
        if self.model not in _MODEL_CODES:
            raise ValueError(f"model must be one of {list(_MODEL_CODES)}")
        if self.redistribution not in _REDIST_CODES:
            raise ValueError("redistribution must be None, 'T' or 'Q'")
        if not (0 <= self.recovery <= 1 and 0 <= self.bankruptcy_cost < 1
                and self.exclusion_period >= 0):
            raise ValueError("need 0<=recovery<=1, 0<=bankruptcy_cost<1, exclusion_period>=0")
        if (self.recovery != 1 or self.bankruptcy_cost != 0 or self.exclusion_period != 0) \
                and (self.model != "L" or self.insolvency != "bankrupt"):
            raise ValueError("bankruptcy options need model='L' and insolvency='bankrupt'")
        if self.insolvency not in _INSOLV_CODES:
            raise ValueError(f"insolvency must be one of {list(_INSOLV_CODES)}")
        if not 0 <= self.savings < 1:
            raise ValueError("savings rate must be in [0, 1)")


@njit(cache=True)
def gini(x):
    """Eq (7): g = 2*sum(i*r_i)/(N*sum(r_i)) - (N+1)/N, r = sorted wealth."""
    r = np.sort(x)
    n = r.shape[0]
    weighted = 0.0
    total = 0.0
    for i in range(n):
        weighted += (i + 1) * r[i]
        total += r[i]
    return 2.0 * weighted / (n * total) - (n + 1.0) / n


@njit(cache=True)
def _transfer_T(m, xi):
    """Eq (4), all agents simultaneously: m_i <- (1-xi) m_i + xi * mean_{j!=i} m_j."""
    n = m.shape[0]
    s = m.sum()
    for i in range(n):
        m[i] = (1.0 - xi) * m[i] + xi * (s - m[i]) / (n - 1)


@njit(cache=True)
def _transfer_Q(m, xi, q):
    """Eq (5)-(6). Threshold is the fraction 1/q of the *maximum wealth*
    (not a population percentile): above it agents pay xi*m_i, at or below it
    agents split the pot equally."""
    n = m.shape[0]
    thr = m.max() / q
    n_q = 0
    s_q = 0.0
    for i in range(n):
        if m[i] > thr:
            s_q += m[i]
        else:
            n_q += 1
    if n_q == 0:
        return
    gift = xi * s_q / n_q
    for i in range(n):
        if m[i] > thr:
            m[i] = (1.0 - xi) * m[i]
        else:
            m[i] = m[i] + gift


@njit(cache=True)
def _run(model, redist, n, lam, rho, dw, xi, tp, q, shared_delta, insolv,
         recovery, cost, excl, sample_times, seed):
    np.random.seed(seed)
    m = np.ones(n)
    n_s = sample_times.shape[0]
    ginis = np.empty(n_s)
    totals = np.empty(n_s)
    minima = np.empty(n_s)
    snaps = np.empty((n_s, n))
    defaults = np.empty(n_s, dtype=np.int64)
    claims = np.empty(n_s)
    paid_at = np.empty(n_s)
    n_def = 0
    claim_sum = 0.0
    paid_sum = 0.0
    blocked = np.zeros(n, dtype=np.int64)  # borrower-ban expiry time per agent
    k = 0
    t_end = sample_times[n_s - 1]
    for t in range(1, t_end + 1):
        i = np.random.randint(0, n)
        if excl > 0 and model == 1:
            while blocked[i] > t:
                i = np.random.randint(0, n)
        j = np.random.randint(0, n - 1)
        if j >= i:
            j += 1  # j != i, uniform over the other agents
        mi = m[i]
        mj = m[j]
        if model == 0:  # R, Eq (1)
            eps = np.random.random()
            pot = (1.0 - lam) * (mi + mj)
            m[i] = lam * mi + eps * pot
            m[j] = lam * mj + (1.0 - eps) * pot
        elif model == 1:  # L, Eq (2): i borrower, j lender
            d = np.random.uniform(-dw, dw)
            new_i = lam * mi + (1.0 - lam) * (mi - rho * mj + d * (mi + mj))
            new_j = lam * mj + (1.0 - lam) * (1.0 + rho) * mj
            if insolv == 3:
                # borrower pays interest only out of what it has left after
                # profit/loss; the lender receives what was actually paid.
                before = lam * mi + (1.0 - lam) * (mi + d * (mi + mj))
                due = (1.0 - lam) * rho * mj
                paid = min(due, max(before, 0.0))
                new_i = max(before - paid, 0.0)
                new_j = mj + paid
            if new_i < 0.0:
                n_def += 1
            if new_i < 0.0 and insolv == 2:
                pass  # trade rejected, both agents keep their wealth
            elif new_i < 0.0 and insolv == 4:
                # Default: limited liability for the borrower, the lender
                # writes off the negative balance against its own wealth.
                claim = new_j - lam * mj
                x = (new_i + claim) * (1.0 - cost)   # assets left to liquidate
                pay = min(recovery * claim, x)
                m[i] = x - pay
                m[j] = lam * mj + pay
                claim_sum += claim
                paid_sum += pay
                if excl > 0:
                    blocked[i] = t + excl
            else:
                if new_i < 0.0 and insolv == 1:
                    new_i = 0.0
                m[i] = new_i
                m[j] = new_j
        else:  # J, Eq (3)
            di = np.random.uniform(-dw, dw)
            dj = di if shared_delta else np.random.uniform(-dw, dw)
            m[i] = lam * mi + (1.0 - lam) * (1.0 + di) * mi
            m[j] = lam * mj + (1.0 - lam) * (1.0 + dj) * mj

        # Record *before* any transfer due at this t: Fig 5's curves meet the
        # exchange-only curves of Fig 4 at t = tp, i.e. g is read just before
        # the transfer fires.
        if k < n_s and t == sample_times[k]:
            ginis[k] = gini(m)
            totals[k] = m.sum()
            minima[k] = m.min()
            defaults[k] = n_def
            claims[k] = claim_sum
            paid_at[k] = paid_sum
            snaps[k, :] = m
            k += 1
        if redist != 0 and t % tp == 0:
            if redist == 1:
                _transfer_T(m, xi)
            else:
                _transfer_Q(m, xi, q)
    return ginis, totals, minima, snaps, defaults, claims, paid_at


@dataclass
class Result:
    times: np.ndarray
    gini: np.ndarray
    total_wealth: np.ndarray
    min_wealth: np.ndarray
    wealth: np.ndarray  # (len(times), N) snapshots
    defaults: np.ndarray  # cumulative borrower defaults (L model; 0 if none/other)
    recovery_rate: np.ndarray  # cumulative lender receipts / claims over defaults (nan if none)


def simulate(p: Params, seed: int = 0) -> Result:
    """Run one realisation and record state at p.samples (sorted, ints >= 1)."""
    times = np.array(sorted(p.samples), dtype=np.int64)
    if times.size == 0 or times[0] < 1:
        raise ValueError("Params.samples must hold at least one time >= 1")
    g, tot, mn, snaps, dfl, clm, pd_ = _run(
        _MODEL_CODES[p.model], _REDIST_CODES[p.redistribution], p.n_agents,
        p.savings, p.interest, p.delta_width, p.transfer_rate,
        p.transfer_period, p.quantile, p.shared_delta, _INSOLV_CODES[p.insolvency],
        p.recovery, p.bankruptcy_cost, p.exclusion_period, times, seed)
    with np.errstate(invalid="ignore", divide="ignore"):
        rec = np.where(clm > 0, pd_ / clm, np.nan)
    return Result(times, g, tot, mn, snaps, dfl, rec)
