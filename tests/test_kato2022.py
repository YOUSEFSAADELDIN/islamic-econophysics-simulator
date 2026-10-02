"""Checks of the Kato (2022) engine: invariants, then the paper's headline results.

Statistical tests use a handful of seeds and generous bands; the full
comparison against the figures lives in experiments/replicate_kato2022.py.
"""

import numpy as np
import pytest

from islamic_econophysics.kato2022 import (
    Params, gini, simulate, _transfer_Q, _transfer_T)

SEEDS = range(6)


def g_at(t, seeds=SEEDS, **kw):
    return float(np.mean([simulate(Params(samples=(t,), **kw), seed=s).gini[0]
                          for s in seeds]))


# ---- Gini (Eq 7) ---------------------------------------------------------
def test_gini_extremes():
    assert gini(np.ones(100)) == pytest.approx(0.0, abs=1e-12)
    x = np.zeros(100); x[0] = 1.0
    assert gini(x) == pytest.approx(1 - 1 / 100)  # (N-1)/N: finite-N delta distribution


def test_gini_exponential_is_half():
    x = np.random.default_rng(0).exponential(size=200_000)
    assert gini(x) == pytest.approx(0.5, abs=0.01)


# ---- invariants ----------------------------------------------------------
def test_r_model_conserves_total_wealth():
    r = simulate(Params(model="R", samples=(10**4,)), seed=1)
    assert r.total_wealth[0] == pytest.approx(1000.0, rel=1e-9)


def test_transfer_T_conserves_wealth_and_contracts():
    m = np.random.default_rng(0).exponential(size=500)
    before = m.sum()
    _transfer_T(m, 0.5)
    assert m.sum() == pytest.approx(before)
    assert gini(m) < gini(np.random.default_rng(0).exponential(size=500))


def test_transfer_Q_conserves_wealth_and_q1_is_noop():
    m = np.random.default_rng(1).exponential(size=500)
    before = m.copy()
    _transfer_Q(m, 0.5, 1.0)  # threshold = m_max: nobody is above it
    assert np.array_equal(m, before)
    _transfer_Q(m, 0.5, 6.0)
    assert m.sum() == pytest.approx(before.sum())
    assert gini(m) < gini(before)


def test_reproducible_by_seed():
    a = simulate(Params(model="J", samples=(5000,)), seed=3).wealth
    b = simulate(Params(model="J", samples=(5000,)), seed=3).wealth
    assert np.array_equal(a, b)


def test_l_model_literal_equations_leave_nonnegative_orthant():
    # Eq (2) taken literally lets the borrower go negative; the default
    # insolvency rule must prevent that.
    lit = simulate(Params(model="L", insolvency="none", samples=(10**5,)), seed=0)
    safe = simulate(Params(model="L", samples=(10**5,)), seed=0)
    assert lit.min_wealth[0] < 0
    assert safe.min_wealth[0] >= 0


def test_bankruptcy_keeps_wealth_nonnegative_and_conserves_total():
    # With no profit/loss noise Eq (2) conserves wealth; defaults happen anyway
    # (rho on a rich lender's wealth can exceed a poor borrower's) and the lender
    # absorbing them must keep the total exact.
    r = simulate(Params(model="L", delta_width=0.0, interest=0.3,
                        samples=(2 * 10**4,)), seed=2)
    assert r.defaults[0] > 0
    assert r.min_wealth[0] >= 0
    assert r.total_wealth[0] == pytest.approx(1000.0, rel=1e-9)


def test_bankruptcy_does_not_inflate_total_wealth():
    # Contrast with "skip", which inflates it by ~2x by t=1e5.
    tot = np.mean([simulate(Params(model="L", samples=(10**5,)), seed=s).total_wealth[0]
                   for s in range(4)])
    assert tot == pytest.approx(1000.0, rel=0.15)


def test_default_counter_zero_without_defaults():
    r = simulate(Params(model="J", samples=(1000,)), seed=0)
    assert r.defaults[0] == 0


# ---- Kato's results (Section "Results") ---------------------------------
def test_R_gini_is_gamma_k2_value():
    # lambda=0.25 -> Chatterjee gamma shape k = 1 + 3*lambda/(1-lambda) = 2,
    # whose Gini is Gamma(k+1/2)/(sqrt(pi) Gamma(k+1)) = 0.375 (Kato: "~0.4").
    assert g_at(10**5, model="R") == pytest.approx(0.375, abs=0.02)


@pytest.mark.parametrize("dw,band", [(0.1, (0.28, 0.40)), (0.2, (0.55, 0.67))])
def test_J_gini_at_1e5(dw, band):
    assert band[0] < g_at(10**5, model="J", delta_width=dw) < band[1]  # Fig 4


def test_J_larger_loss_width_is_more_unequal():  # Kato: gharar ban
    assert g_at(10**5, model="J", delta_width=0.2) > g_at(10**5, model="J", delta_width=0.1) + 0.15


def test_J_heads_to_one():
    assert g_at(10**6, seeds=range(3), model="J", delta_width=0.2) > 0.95


def test_L_more_unequal_than_J_same_width():  # Kato: one party bears all risk
    assert g_at(10**5, model="L", interest=0.0) > g_at(10**5, model="J", delta_width=0.1) + 0.1


def test_interest_free_L_still_concentrates():  # the "risk, not interest" point
    assert g_at(10**5, model="L", interest=0.0) > 0.4


@pytest.mark.parametrize("model,kw,tp,band", [
    ("J", dict(delta_width=0.1), 10**5, (0.33, 0.43)),
    ("J", dict(delta_width=0.2), 10**5, (0.60, 0.72)),
    ("J", dict(delta_width=0.1), 10**4, (0.09, 0.17)),
    ("J", dict(delta_width=0.2), 10**4, (0.20, 0.30)),
    ("L", dict(interest=0.05), 10**4, (0.26, 0.37)),
])
def test_T_model_plateaus_fig5(model, kw, tp, band):
    g = g_at(10**6, model=model, redistribution="T", transfer_rate=0.5,
             transfer_period=tp, **kw)
    assert band[0] < g < band[1]


def test_shorter_transfer_period_lowers_gini():
    kw = dict(model="J", delta_width=0.2, redistribution="T", transfer_rate=0.5)
    assert g_at(10**6, transfer_period=10**4, **kw) < g_at(10**6, transfer_period=10**5, **kw) - 0.2


def test_gini_falls_with_xi_then_saturates():  # Fig 6
    kw = dict(model="J", delta_width=0.1, redistribution="T", transfer_period=10**4)
    g0, g2, g5, g10 = (g_at(10**6, transfer_rate=x, **kw) for x in (0.0, 0.2, 0.5, 1.0))
    assert g0 - g2 > 0.5          # steep fall up to xi ~ 0.2
    assert abs(g2 - g5) < 0.06    # then it flattens
    assert abs(g5 - g10) < 0.04


def test_Q_matches_T_at_mid_quantile_for_J01():  # Fig 7, "1/4 to 1/6 quantile"
    kw = dict(model="J", delta_width=0.1, transfer_rate=0.5, transfer_period=10**5)
    gt = g_at(10**6, redistribution="T", **kw)
    for q in (5.0, 6.0):
        assert g_at(10**6, redistribution="Q", quantile=q, **kw) == pytest.approx(gt, abs=0.05)


def test_R_T_is_almost_R():
    # Kato: "R-T results are almost identical to R, regardless of redistribution".
    r = g_at(10**5, model="R")
    rt = g_at(10**5, model="R", redistribution="T", transfer_period=10**4)
    assert rt == pytest.approx(r, abs=0.03)
