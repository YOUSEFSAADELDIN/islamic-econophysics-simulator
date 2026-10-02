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
    #   lender_exclusion_period E: ... and cannot lend for the next E steps
    #   garnishment_rate g > 0: debt ledger. The unpaid part U of the claim becomes a debt
    #                      of the borrower to that lender (one open debt per borrower; a
    #                      later lender of an already indebted borrower just absorbs its loss).
    #                      Whenever the debtor's wealth rises (profit as borrower, interest as
    #                      lender, transfers received) a fraction g of the gain goes to the
    #                      creditor until U is repaid. Needs recovery=1 and bankruptcy_cost=0.
    #                      g=1 is the closest to Eq (2), where the shortfall stays as negative wealth.
    recovery: float = 1.0
    bankruptcy_cost: float = 0.0
    exclusion_period: int = 0
    lender_exclusion_period: int = 0
    garnishment_rate: float = 0.0
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
                and self.exclusion_period >= 0 and self.lender_exclusion_period >= 0
                and 0 <= self.garnishment_rate <= 1):
            raise ValueError("need 0<=recovery<=1, 0<=bankruptcy_cost<1, exclusion periods>=0, "
                             "0<=garnishment_rate<=1")
        opts = (self.recovery != 1 or self.bankruptcy_cost != 0 or self.exclusion_period != 0
                or self.lender_exclusion_period != 0 or self.garnishment_rate != 0)
        if opts and (self.model != "L" or self.insolvency != "bankrupt"):
            raise ValueError("bankruptcy options need model='L' and insolvency='bankrupt'")
        if self.garnishment_rate > 0 and (self.recovery != 1 or self.bankruptcy_cost != 0):
            raise ValueError("the debt ledger replaces recovery/bankruptcy_cost; use them separately")
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
def _garnish_all(m, before, debt_amt, debt_cred, g):
    """Divert a fraction g of every indebted agent's wealth *gain* since `before`
    to its creditor (gains are measured first, then paid). Returns total repaid."""
    n = m.shape[0]
    gains = np.zeros(n)
    for a in range(n):
        if debt_amt[a] > 0.0 and m[a] > before[a]:
            gains[a] = m[a] - before[a]
    repaid = 0.0
    for a in range(n):
        if gains[a] > 0.0:
            pay = min(g * gains[a], debt_amt[a])
            m[a] -= pay
            m[debt_cred[a]] += pay
            debt_amt[a] -= pay
            repaid += pay
    return repaid


_MAX_TRIES = 10000  # bounded resampling when many agents are banned


@njit(cache=True)
def _run(model, redist, n, lam, rho, dw, xi, tp, q, shared_delta, insolv,
         recovery, cost, excl, lend_excl, garnish, sample_times, seed):
    np.random.seed(seed)
    m = np.ones(n)
    n_s = sample_times.shape[0]
    ginis = np.empty(n_s)
    totals = np.empty(n_s)
    minima = np.empty(n_s)
    snaps = np.empty((n_s, n))
    # per-sample diagnostics: defaults, claims, paid at default, debt outstanding,
    # debt repaid, idle steps (no eligible pair)
    diag = np.zeros((n_s, 6))
    n_def = 0
    claim_sum = 0.0
    paid_sum = 0.0
    repaid_sum = 0.0
    idle = 0
    blocked = np.zeros(n, dtype=np.int64)       # borrower-ban expiry time per agent
    lend_blocked = np.zeros(n, dtype=np.int64)  # lender-ban expiry time per agent
    debt_amt = np.zeros(n)                      # open debt owed by agent
    debt_cred = np.full(n, -1, dtype=np.int64)  # ... and to whom
    k = 0
    t_end = sample_times[n_s - 1]
    for t in range(1, t_end + 1):
        ok = True
        i = np.random.randint(0, n)
        if excl > 0 and model == 1:
            tries = 0
            while blocked[i] > t and tries < _MAX_TRIES:
                i = np.random.randint(0, n)
                tries += 1
            if blocked[i] > t:
                ok = False
        j = np.random.randint(0, n - 1)
        if j >= i:
            j += 1  # j != i, uniform over the other agents
        if lend_excl > 0 and model == 1 and ok:
            tries = 0
            while lend_blocked[j] > t and tries < _MAX_TRIES:
                j = np.random.randint(0, n - 1)
                if j >= i:
                    j += 1
                tries += 1
            if lend_blocked[j] > t:
                ok = False
        if not ok:
            idle += 1
        else:
            mi = m[i]
            mj = m[j]
            defaulted = False
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
                    defaulted = True
                    if excl > 0:
                        blocked[i] = t + excl
                    if lend_excl > 0:
                        lend_blocked[i] = t + lend_excl
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

            if garnish > 0.0:
                # Income garnishment for indebted members of the pair.
                gi = m[i] - mi
                gj = m[j] - mj
                if debt_amt[i] > 0.0 and gi > 0.0:
                    pay = min(garnish * gi, debt_amt[i])
                    m[i] -= pay
                    m[debt_cred[i]] += pay
                    debt_amt[i] -= pay
                    repaid_sum += pay
                if debt_amt[j] > 0.0 and gj > 0.0:
                    pay = min(garnish * gj, debt_amt[j])
                    m[j] -= pay
                    m[debt_cred[j]] += pay
                    debt_amt[j] -= pay
                    repaid_sum += pay
                if defaulted and debt_amt[i] <= 0.0:
                    # unpaid part of the claim becomes a debt to this lender
                    debt_amt[i] = -new_i
                    debt_cred[i] = j

        # Record *before* any transfer due at this t: Fig 5's curves meet the
        # exchange-only curves of Fig 4 at t = tp, i.e. g is read just before
        # the transfer fires.
        if k < n_s and t == sample_times[k]:
            ginis[k] = gini(m)
            totals[k] = m.sum()
            minima[k] = m.min()
            diag[k, 0] = n_def
            diag[k, 1] = claim_sum
            diag[k, 2] = paid_sum
            diag[k, 3] = debt_amt.sum()
            diag[k, 4] = repaid_sum
            diag[k, 5] = idle
            snaps[k, :] = m
            k += 1
        if redist != 0 and t % tp == 0:
            if garnish > 0.0:
                before_t = m.copy()
            if redist == 1:
                _transfer_T(m, xi)
            else:
                _transfer_Q(m, xi, q)
            if garnish > 0.0:
                repaid_sum += _garnish_all(m, before_t, debt_amt, debt_cred, garnish)
    return ginis, totals, minima, snaps, diag


@dataclass
class Result:
    times: np.ndarray
    gini: np.ndarray
    total_wealth: np.ndarray
    min_wealth: np.ndarray
    wealth: np.ndarray  # (len(times), N) snapshots
    defaults: np.ndarray  # cumulative borrower defaults (L model; 0 if none/other)
    recovery_rate: np.ndarray  # lender receipts (at default + later repayments) / claims over defaults (nan if none)
    debt_outstanding: np.ndarray  # open ledger debt (garnishment_rate > 0)
    debt_repaid: np.ndarray    # cumulative ledger repayments
    idle_steps: np.ndarray     # steps with no eligible pair (everyone banned)


def simulate(p: Params, seed: int = 0) -> Result:
    """Run one realisation and record state at p.samples (sorted, ints >= 1)."""
    times = np.array(sorted(p.samples), dtype=np.int64)
    if times.size == 0 or times[0] < 1:
        raise ValueError("Params.samples must hold at least one time >= 1")
    g, tot, mn, snaps, diag = _run(
        _MODEL_CODES[p.model], _REDIST_CODES[p.redistribution], p.n_agents,
        p.savings, p.interest, p.delta_width, p.transfer_rate,
        p.transfer_period, p.quantile, p.shared_delta, _INSOLV_CODES[p.insolvency],
        p.recovery, p.bankruptcy_cost, p.exclusion_period, p.lender_exclusion_period,
        p.garnishment_rate, times, seed)
    defaults, claims, paid, debt_out, repaid, idle = diag.T
    with np.errstate(invalid="ignore", divide="ignore"):
        rec = np.where(claims > 0, (paid + repaid) / claims, np.nan)
    return Result(times, g, tot, mn, snaps, defaults.astype(np.int64), rec,
                  debt_out, repaid, idle.astype(np.int64))
