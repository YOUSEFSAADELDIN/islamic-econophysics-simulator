"""Yard-sale model with wealth-attained advantage (WAA) and a zakat levy.

Chakraborti (2002) / Hayes yard-sale model, extended with the WAA of Boghosian et al.
(Physica A 476, 2017) and a periodic zakat (wealth levy with a nisab threshold). This is
independent of Kato (2022): it asks whether a zakat-style levy stops the oligarchy that a
wealth-biased yard sale produces.

Trade (one random pair per step):
    stake  = beta * min(w_i, w_j)                        (the poorer agent's wealth is at risk)
    p_i    = 1/2 * (1 + waa * (w_i - w_j) / (w_i + w_j)) (probability that i wins)
    the winner gains `stake`, the loser pays it. waa = 0 is the fair yard sale.
Zakat (every `trades_per_agent_year * N / 2` trades = one lunar year):
    every agent with w >= nisab * mean(w) pays zakat_rate * w
    the pool goes to all agents equally ("all") or to the agents below the nisab ("poor").

Recipients "eight": the pool is split between the eight Qur'anic categories (Q 9:60) by
`category_weights`. A yard-sale economy has no slaves, travellers or loan book, so categories
4-8 are *modelling choices*, listed in CATEGORIES and docs/zakat_categories.md.

The mapping from the continuous-time Fokker-Planck parameters (chi, zeta) of Boghosian et al.
to this discrete model is *not* assumed: `waa` is the bias of the coin, and the critical zakat
rate is measured by simulation (experiments/zakat_sweep.py).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from .kato2022 import gini

RECIPIENTS = {"all": 0, "poor": 1, "eight": 2}

CATEGORIES = (
    "fuqara",      # the poor: wealth < poverty_line x mean; paid pro rata to the gap to the nisab
    "masakin",     # the needy: poverty_line <= wealth < nisab; same rule
    "amilin",      # administrators: wages, equal split over collector_frac x N random agents
    "muallafat",   # [interpretive] outreach: equal split over muallaf_frac x N random agents below mean
    "riqab",       # [interpretive] release from persistent poverty: agents below the poverty line
                   # for >= riqab_years consecutive years, equal split
    "gharimin",    # [interpretive] distressed debtors: lost > gharimin_loss of last year's wealth and
                   # now below the nisab; paid pro rata to the loss (never more than the loss)
    "sabilillah",  # public good: equal per-capita grant to everyone (also takes any unspent share)
    "ibn_sabil",   # [interpretive] stranded wayfarer: wayfarer_frac x N random agents of any wealth
)

# Stylised allocation profiles (weights over CATEGORIES). They are NOT a fiqh ruling: see the doc.
PROFILES = {
    "all_public":      (0, 0, 0, 0, 0, 0, 1, 0),          # = the earlier recipients="all"
    "poor_only":       (0.5, 0.5, 0, 0, 0, 0, 0, 0),
    "eight_equal":     (1 / 8,) * 8,                      # every category present gets 1/8
    "seven_equal":     (1 / 7, 1 / 7, 1 / 7, 0, 1 / 7, 1 / 7, 1 / 7, 1 / 7),  # no muallafat
    "needy_first":     (0.35, 0.35, 0.1, 0, 0.05, 0.1, 0.05, 0),
}


@dataclass(frozen=True)
class YSParams:
    n_agents: int = 500
    beta: float = 0.1                    # fraction of the poorer agent's wealth staked per trade
    waa: float = 0.0                     # wealth-attained advantage, 0 <= waa <= 1
    zakat_rate: float = 0.0              # per lunar year, on w for agents above the nisab
    nisab: float = 0.0                   # threshold as a fraction of mean wealth
    recipients: str = "all"              # "all" or "poor" (agents below the nisab)
    trades_per_agent_year: float = 20.0  # how many trades each agent makes per lunar year
    years: int = 1000
    n_samples: int = 100                 # equally spaced in time, last one at `years`
    # recipients="eight" only (needs nisab > 0):
    category_weights: tuple = PROFILES["eight_equal"]
    poverty_line: float = 0.2            # x mean wealth; also the poverty line of the headcount metric
    riqab_years: int = 5
    gharimin_loss: float = 0.5
    collector_frac: float = 0.01
    muallaf_frac: float = 0.05
    wayfarer_frac: float = 0.02

    def __post_init__(self):
        if not 0 <= self.waa <= 1:
            raise ValueError("waa must be in [0, 1]")
        if not 0 < self.beta < 1:
            raise ValueError("beta must be in (0, 1)")
        if not 0 <= self.zakat_rate < 1 or self.nisab < 0:
            raise ValueError("need 0 <= zakat_rate < 1 and nisab >= 0")
        if self.recipients not in RECIPIENTS:
            raise ValueError(f"recipients must be one of {list(RECIPIENTS)}")
        if self.years < 1 or self.n_samples < 1 or self.n_agents < 3:
            raise ValueError("years, n_samples >= 1 and n_agents >= 3")
        if len(self.category_weights) != 8 or min(self.category_weights) < 0 or sum(self.category_weights) <= 0:
            raise ValueError("category_weights: 8 non-negative weights with a positive sum")
        if self.recipients == "eight" and self.nisab <= 0:
            raise ValueError("recipients='eight' needs nisab > 0 (it defines who is needy)")
        if not 0 <= self.poverty_line < max(self.nisab, 1e-12) and self.recipients == "eight":
            raise ValueError("poverty_line must lie below the nisab")
        if self.year_trades < 1:
            raise ValueError("trades_per_agent_year * n_agents / 2 must be >= 1")

    @property
    def year_trades(self) -> int:
        return int(round(self.trades_per_agent_year * self.n_agents / 2))


@njit(cache=True)
def _zakat(w, z, nisab, recip):
    n = w.shape[0]
    thr = nisab * w.sum() / n
    pool = 0.0
    for a in range(n):
        if w[a] >= thr:
            t = z * w[a]
            w[a] -= t
            pool += t
    if recip == 1:
        k = 0
        for a in range(n):
            if w[a] < thr:
                k += 1
        if k > 0:
            for a in range(n):
                if w[a] < thr:
                    w[a] += pool / k
            return pool
    for a in range(n):
        w[a] += pool / n
    return pool


@njit(cache=True)
def _by_need(w, got, mask, need, total):
    """Pay `total` to the agents in `mask` pro rata to `need`, never more than `need`.
    Returns the unspent part."""
    s = 0.0
    for a in range(w.shape[0]):
        if mask[a] and need[a] > 0.0:
            s += need[a]
    if s <= 0.0:
        return total
    paid = min(total, s)
    for a in range(w.shape[0]):
        if mask[a] and need[a] > 0.0:
            p = paid * need[a] / s
            w[a] += p
            got[a] += p
    return total - paid


@njit(cache=True)
def _equal_to(w, got, idx, k, total):
    for q in range(k):
        w[idx[q]] += total / k
        got[idx[q]] += total / k


@njit(cache=True)
def _random_pick(mask, k):
    """Indices of up to k agents drawn uniformly without replacement from `mask`."""
    n = mask.shape[0]
    keys = np.random.random(n)
    cnt = 0
    for a in range(n):
        if mask[a]:
            cnt += 1
        else:
            keys[a] = 2.0
    kk = min(k, cnt)
    return np.argsort(keys)[:kk], kk


@njit(cache=True)
def _zakat_eight(w, w_prev, pov_years, z, nisab, weights, cat, rcv):
    """Levy as in _zakat, then split the pool over the eight categories. cat = [poverty_line,
    riqab_years, gharimin_loss, collector_frac, muallaf_frac, wayfarer_frac]. rcv[k] accumulates
    what category k actually paid out; rcv[8] the total pool, rcv[9] what the poorest 20% of
    agents received (all categories, public share included) and rcv[10] the unspent part of
    categories with nobody eligible, which is returned to everyone equally. Returns the pool."""
    n = w.shape[0]
    mean = w.sum() / n
    thr = nisab * mean
    pool = 0.0
    for a in range(n):
        if w[a] >= thr:
            t = z * w[a]
            w[a] -= t
            pool += t
    pl = cat[0] * mean
    w0 = w.copy()      # who is eligible for what is decided on wealth before any grant
    got = np.zeros(n)
    order = np.argsort(w)
    bot = np.zeros(n, dtype=np.bool_)
    for r in range(n // 5):
        bot[order[r]] = True
    public = 0.0
    unspent = 0.0
    wsum = weights.sum()
    mask = np.zeros(n, dtype=np.bool_)
    need = np.zeros(n)
    allmask = np.ones(n, dtype=np.bool_)
    for k in range(8):
        amt = pool * weights[k] / wsum
        if amt <= 0.0:
            continue
        left = 0.0
        if k == 0 or k == 1:
            for a in range(n):
                mask[a] = (w0[a] < pl) if k == 0 else (w0[a] >= pl and w0[a] < thr)
                need[a] = thr - w0[a]
            left = _by_need(w, got, mask, need, amt)
        elif k == 2:
            idx, kk = _random_pick(allmask, max(1, int(cat[3] * n)))
            _equal_to(w, got, idx, kk, amt)
        elif k == 3:
            for a in range(n):
                mask[a] = w0[a] < mean
            idx, kk = _random_pick(mask, max(1, int(cat[4] * n)))
            if kk == 0:
                left = amt
            else:
                _equal_to(w, got, idx, kk, amt)
        elif k == 4:
            kk = 0
            for a in range(n):
                mask[a] = pov_years[a] >= cat[1]
                if mask[a]:
                    kk += 1
            if kk == 0:
                left = amt
            else:
                for a in range(n):
                    if mask[a]:
                        w[a] += amt / kk
                        got[a] += amt / kk
        elif k == 5:
            for a in range(n):
                loss = w_prev[a] - w0[a]
                mask[a] = w_prev[a] > 0.0 and loss > cat[2] * w_prev[a] and w0[a] < thr
                need[a] = loss
            left = _by_need(w, got, mask, need, amt)
        elif k == 6:
            public += amt
        else:
            idx, kk = _random_pick(allmask, max(1, int(cat[5] * n)))
            _equal_to(w, got, idx, kk, amt)
        rcv[k] += amt - left
        public += left
        unspent += left
    if public > 0.0:
        for a in range(n):
            w[a] += public / n
            got[a] += public / n
    rcv[10] += unspent
    rcv[8] += pool
    for a in range(n):
        if bot[a]:
            rcv[9] += got[a]
    # trackers for next year
    for a in range(n):
        if w[a] < pl:
            pov_years[a] += 1
        else:
            pov_years[a] = 0
        w_prev[a] = w[a]
    return pool


@njit(cache=True)
def _run(n, beta, waa, z, nisab, recip, year_trades, years, n_samples, seed, weights, cat):
    np.random.seed(seed)
    w = np.ones(n)
    w_prev = np.ones(n)
    pov_years = np.zeros(n, dtype=np.int64)
    rcv = np.zeros(11)
    g = np.empty(n_samples)
    top1 = np.empty(n_samples)
    top10 = np.empty(n_samples)
    over = np.empty(n_samples)
    pov = np.empty(n_samples)
    pov_pre = np.empty(n_samples)
    bot50 = np.empty(n_samples)
    sample_years = np.empty(n_samples, dtype=np.int64)
    for s in range(n_samples):
        sample_years[s] = max(1, (years * (s + 1)) // n_samples)
    s = 0
    for y in range(1, years + 1):
        for _ in range(year_trades):
            i = np.random.randint(0, n)
            j = np.random.randint(0, n - 1)
            if j >= i:
                j += 1
            wi = w[i]
            wj = w[j]
            stake = beta * min(wi, wj)
            p = 0.5 * (1.0 + waa * (wi - wj) / (wi + wj))
            if np.random.random() < p:
                w[i] = wi + stake
                w[j] = wj - stake
            else:
                w[i] = wi - stake
                w[j] = wj + stake
        if s < n_samples and y == sample_years[s]:
            # poverty at the end of the trading year, *before* zakat arrives (worst point of the cycle)
            tot0 = w.sum()
            c = 0
            for a in range(n):
                if w[a] < cat[0] * tot0 / n:
                    c += 1
            pov_pre[s] = c / n
        if z > 0.0:
            if recip == 2:
                pool = _zakat_eight(w, w_prev, pov_years, z, nisab, weights, cat, rcv)
            else:
                pool = _zakat(w, z, nisab, recip)
                rcv[8] += pool
        while s < n_samples and y == sample_years[s]:
            tot = w.sum()
            srt = np.sort(w)
            g[s] = gini(w)
            top1[s] = srt[n - 1] / tot
            k = max(1, n // 10)
            top10[s] = srt[n - k:].sum() / tot
            bot50[s] = srt[: n // 2].sum() / tot
            thr = nisab * tot / n
            c = 0
            for a in range(n):
                if w[a] >= thr:
                    c += 1
            over[s] = c / n
            c = 0
            for a in range(n):
                if w[a] < cat[0] * tot / n:
                    c += 1
            pov[s] = c / n
            s += 1
    return sample_years, g, top1, top10, over, w, pov, bot50, rcv, pov_pre


@dataclass
class YSResult:
    years: np.ndarray
    gini: np.ndarray
    top1: np.ndarray        # wealth share of the richest agent (the "oligarch")
    top10pct: np.ndarray    # wealth share of the richest 10% of agents
    frac_over_nisab: np.ndarray
    wealth: np.ndarray      # final wealth vector
    poverty: np.ndarray     # share of agents below poverty_line x mean wealth, just after zakat
    poverty_pre: np.ndarray  # same, at the end of the trading year before zakat (worst point of the cycle)
    bottom50: np.ndarray    # wealth share of the poorer half
    paid: np.ndarray        # cumulative zakat statistics, see `stats`

    def steady(self, field: str, tail: float = 0.25) -> float:
        """Mean of `field` over the last `tail` fraction of the samples."""
        x = getattr(self, field)
        k = max(1, int(len(x) * tail))
        return float(np.mean(x[-k:]))

    def stats(self) -> dict:
        """Cumulative zakat accounting over the whole run (recipients='eight' fills all of it;
        other rules only the total). Shares are fractions of the total pool collected."""
        r = self.paid
        pool = r[8]
        out = {"pool": float(pool)}
        if pool > 0 and r[:8].sum() + r[10] > 0:
            for k, name in enumerate(CATEGORIES):
                out[f"share_{name}"] = float(r[k] / pool)
            out["share_unspent_returned"] = float(r[10] / pool)
            out["share_to_poorest_20pct"] = float(r[9] / pool)
        return out


def simulate(p: YSParams, seed: int = 0) -> YSResult:
    wts = np.array(p.category_weights, dtype=np.float64)
    cat = np.array([p.poverty_line, p.riqab_years, p.gharimin_loss, p.collector_frac,
                    p.muallaf_frac, p.wayfarer_frac], dtype=np.float64)
    yrs, g, t1, t10, ov, w, pov, b50, rcv, pov_pre = _run(
        p.n_agents, p.beta, p.waa, p.zakat_rate, p.nisab, RECIPIENTS[p.recipients],
        p.year_trades, p.years, p.n_samples, seed, wts, cat)
    return YSResult(yrs, g, t1, t10, ov, w, pov, pov_pre, b50, rcv)
