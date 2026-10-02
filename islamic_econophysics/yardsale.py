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

The mapping from the continuous-time Fokker-Planck parameters (chi, zeta) of Boghosian et al.
to this discrete model is *not* assumed: `waa` is the bias of the coin, and the critical zakat
rate is measured by simulation (experiments/zakat_sweep.py).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from .kato2022 import gini

RECIPIENTS = {"all": 0, "poor": 1}


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
            return
    for a in range(n):
        w[a] += pool / n


@njit(cache=True)
def _run(n, beta, waa, z, nisab, recip, year_trades, years, n_samples, seed):
    np.random.seed(seed)
    w = np.ones(n)
    g = np.empty(n_samples)
    top1 = np.empty(n_samples)
    top10 = np.empty(n_samples)
    over = np.empty(n_samples)
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
        if z > 0.0:
            _zakat(w, z, nisab, recip)
        while s < n_samples and y == sample_years[s]:
            tot = w.sum()
            srt = np.sort(w)
            g[s] = gini(w)
            top1[s] = srt[n - 1] / tot
            k = max(1, n // 10)
            top10[s] = srt[n - k:].sum() / tot
            thr = nisab * tot / n
            c = 0
            for a in range(n):
                if w[a] >= thr:
                    c += 1
            over[s] = c / n
            s += 1
    return sample_years, g, top1, top10, over, w


@dataclass
class YSResult:
    years: np.ndarray
    gini: np.ndarray
    top1: np.ndarray        # wealth share of the richest agent (the "oligarch")
    top10pct: np.ndarray    # wealth share of the richest 10% of agents
    frac_over_nisab: np.ndarray
    wealth: np.ndarray      # final wealth vector

    def steady(self, field: str, tail: float = 0.25) -> float:
        """Mean of `field` over the last `tail` fraction of the samples."""
        x = getattr(self, field)
        k = max(1, int(len(x) * tail))
        return float(np.mean(x[-k:]))


def simulate(p: YSParams, seed: int = 0) -> YSResult:
    yrs, g, t1, t10, ov, w = _run(p.n_agents, p.beta, p.waa, p.zakat_rate, p.nisab,
                                  RECIPIENTS[p.recipients], p.year_trades, p.years,
                                  p.n_samples, seed)
    return YSResult(yrs, g, t1, t10, ov, w)
