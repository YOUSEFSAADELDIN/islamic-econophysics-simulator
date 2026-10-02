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
    #   ledger extensions (need garnishment_rate > 0):
    #     multi_debt:         a borrower may owe several lenders at once (an N x N matrix of
    #                         debts); repayments are split pro rata over the creditors, and a
    #                         new default always adds to what is owed to that lender
    #     seize_receivables:  the defaulter's estate includes the claims others owe it: its
    #                         lender takes them over pro rata up to the unpaid amount (set-off
    #                         if the lender itself owes the defaulter). Needs multi_debt.
    #     debt_term T > 0:    a debt is forgiven (lender writes it off) T steps after the most
    #                         recent default on that lender; 0 = never
    recovery: float = 1.0
    bankruptcy_cost: float = 0.0
    exclusion_period: int = 0
    lender_exclusion_period: int = 0
    garnishment_rate: float = 0.0
    multi_debt: bool = False
    seize_receivables: bool = False
    debt_term: int = 0
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
        if self.debt_term < 0:
            raise ValueError("debt_term must be >= 0")
        if (self.multi_debt or self.seize_receivables or self.debt_term) and self.garnishment_rate <= 0:
            raise ValueError("multi_debt / seize_receivables / debt_term extend the ledger: "
                             "set garnishment_rate > 0")
        if self.seize_receivables and not self.multi_debt:
            raise ValueError("seize_receivables needs multi_debt=True")
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


_EPS = 1e-12
MATERIAL = 0.01  # a default counts as material if the unpaid part exceeds 1% of initial mean wealth


@njit(cache=True)
def _expire_row(a, t, term, D, born, dtot):
    """Forgive agent a's debts older than `term` steps; refresh dtot[a]. Returns amount forgiven."""
    if term <= 0 or dtot[a] <= _EPS:
        return 0.0
    forgiven = 0.0
    remaining = 0.0
    for k in range(D.shape[1]):
        if D[a, k] > 0.0:
            if t - born[a, k] > term:
                forgiven += D[a, k]
                D[a, k] = 0.0
            else:
                remaining += D[a, k]
    dtot[a] = remaining
    return forgiven


@njit(cache=True)
def _garnish_agent(a, gain, g, t, term, m, D, born, dtot):
    """Divert a fraction g of agent a's wealth gain to its creditors, pro rata to what is
    owed to each. Returns (repaid, forgiven) where forgiven = debts that expired first."""
    forgiven = _expire_row(a, t, term, D, born, dtot)
    if gain <= 0.0 or dtot[a] <= _EPS:
        return 0.0, forgiven
    n = D.shape[1]
    tot = 0.0
    for k in range(n):
        tot += D[a, k]
    if tot <= _EPS:
        dtot[a] = 0.0
        return 0.0, forgiven
    pay = min(g * gain, tot)
    frac = pay / tot
    for k in range(n):
        if D[a, k] > 0.0:
            p = D[a, k] * frac
            D[a, k] -= p
            m[k] += p
    m[a] -= pay
    dtot[a] = tot - pay
    return pay, forgiven


@njit(cache=True)
def _garnish_all(m, before, g, t, term, D, born, dtot):
    """Garnish every indebted agent's gain since `before` (gains measured first, then paid).
    Returns (repaid, forgiven)."""
    n = m.shape[0]
    gains = np.zeros(n)
    for a in range(n):
        if m[a] > before[a]:
            gains[a] = m[a] - before[a]
    repaid = 0.0
    forgiven = 0.0
    for a in range(n):
        r, f = _garnish_agent(a, gains[a], g, t, term, m, D, born, dtot)
        repaid += r
        forgiven += f
    return repaid, forgiven


@njit(cache=True)
def _seize_receivables(i, j, u, t, term, D, born, dtot):
    """Defaulter i owes lender j an unpaid amount u. j takes over claims owed to i, pro rata,
    up to u (if j itself owes i, that part is set off). Returns (seized, forgiven)."""
    n = D.shape[0]
    forgiven = 0.0
    rtot = 0.0
    for k in range(n):
        if D[k, i] > 0.0:
            if term > 0 and t - born[k, i] > term:
                forgiven += D[k, i]
                dtot[k] -= D[k, i]
                D[k, i] = 0.0
            else:
                rtot += D[k, i]
    if rtot <= _EPS or u <= 0.0:
        return 0.0, forgiven
    f = min(1.0, u / rtot)
    for k in range(n):
        if D[k, i] > 0.0:
            x = D[k, i] * f
            D[k, i] -= x
            if k == j:
                dtot[j] -= x                    # set-off: j's debt to i cancels
            else:
                D[k, j] += x                    # claim assigned from i to j; k's total unchanged
                born[k, j] = max(born[k, j], born[k, i])
    return f * rtot, forgiven


_MAX_TRIES = 10000  # bounded resampling when many agents are banned


@njit(cache=True)
def _run(model, redist, n, lam, rho, dw, xi, tp, q, shared_delta, insolv,
         recovery, cost, excl, lend_excl, garnish, multi, seize, term,
         sample_times, seed):
    np.random.seed(seed)
    m = np.ones(n)
    n_s = sample_times.shape[0]
    ginis = np.empty(n_s)
    totals = np.empty(n_s)
    minima = np.empty(n_s)
    snaps = np.empty((n_s, n))
    # per-sample diagnostics: defaults, claims, paid at default, debt outstanding,
    # debt repaid, idle steps (no eligible pair), debt forgiven, receivables seized,
    # material defaults (unpaid > MATERIAL), share of agents with negative net worth
    diag = np.zeros((n_s, 10))
    n_def = 0
    claim_sum = 0.0
    paid_sum = 0.0
    repaid_sum = 0.0
    forgiven_sum = 0.0
    seized_sum = 0.0
    n_material = 0
    idle = 0
    blocked = np.zeros(n, dtype=np.int64)       # borrower-ban expiry time per agent
    lend_blocked = np.zeros(n, dtype=np.int64)  # lender-ban expiry time per agent
    nd = n if garnish > 0.0 else 1              # ledger matrices only when the ledger is on
    D = np.zeros((nd, nd))                      # D[a, k]: what agent a owes agent k
    born = np.zeros((nd, nd), dtype=np.int64)   # step of the latest default behind D[a, k]
    dtot = np.zeros(nd)                         # row sums of D (fast "is indebted?" filter)
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
                    if -new_i > MATERIAL:
                        n_material += 1
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
                # Income garnishment for indebted members of the pair (gains measured
                # before either pays, so repayments received are not themselves garnished).
                gi = m[i] - mi
                gj = m[j] - mj
                r1, f1 = _garnish_agent(i, gi, garnish, t, term, m, D, born, dtot)
                r2, f2 = _garnish_agent(j, gj, garnish, t, term, m, D, born, dtot)
                repaid_sum += r1 + r2
                forgiven_sum += f1 + f2
                if defaulted:
                    u = -new_i                     # unpaid part of the claim
                    forgiven_sum += _expire_row(i, t, term, D, born, dtot)
                    if seize:
                        sz, fg = _seize_receivables(i, j, u, t, term, D, born, dtot)
                        u -= sz
                        seized_sum += sz
                        forgiven_sum += fg
                    # Single-debt mode: a borrower with an open debt cannot take on a
                    # second one, so this lender just absorbs the loss.
                    if u > _EPS and (multi or dtot[i] <= _EPS):
                        D[i, j] += u
                        born[i, j] = t
                        dtot[i] += u

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
            if garnish > 0.0 and term > 0:
                for a in range(n):               # make "outstanding" exclude expired debts
                    forgiven_sum += _expire_row(a, t, term, D, born, dtot)
            diag[k, 3] = D.sum() if garnish > 0.0 else 0.0
            diag[k, 4] = repaid_sum
            diag[k, 5] = idle
            diag[k, 6] = forgiven_sum
            diag[k, 7] = seized_sum
            diag[k, 8] = n_material
            if garnish > 0.0:
                # net worth = holdings + claims on others - own debts
                neg = 0
                for a in range(n):
                    recv = 0.0
                    for c in range(n):
                        recv += D[c, a]
                    if m[a] + recv - D[a].sum() < 0.0:
                        neg += 1
                diag[k, 9] = neg / n
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
                r_t, f_t = _garnish_all(m, before_t, garnish, t, term, D, born, dtot)
                repaid_sum += r_t
                forgiven_sum += f_t
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
    debt_forgiven: np.ndarray  # cumulative ledger debt expired under debt_term
    claims_seized: np.ndarray  # cumulative receivables taken over / set off at defaults
    material_defaults: np.ndarray  # defaults with unpaid part > MATERIAL (`defaults` also counts trivial ones between near-zero agents)
    neg_net_worth_share: np.ndarray  # share of agents whose holdings + receivables - debts < 0 (ledger only)


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
        p.garnishment_rate, p.multi_debt, p.seize_receivables, p.debt_term, times, seed)
    defaults, claims, paid, debt_out, repaid, idle, forgiven, seized, material, negnet = diag.T
    with np.errstate(invalid="ignore", divide="ignore"):
        rec = np.where(claims > 0, (paid + repaid) / claims, np.nan)
    return Result(times, g, tot, mn, snaps, defaults.astype(np.int64), rec,
                  debt_out, repaid, idle.astype(np.int64), forgiven, seized,
                  material.astype(np.int64), negnet)
