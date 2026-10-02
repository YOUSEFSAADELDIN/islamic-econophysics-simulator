"""Checks of the Kato (2022) engine: invariants, then the paper's headline results.

Statistical tests use a handful of seeds and generous bands; the full
comparison against the figures lives in experiments/replicate_kato2022.py.
"""

import numpy as np
import pytest

from islamic_econophysics.kato2022 import (
    Params, gini, simulate, _garnish_all, _seize_receivables, _transfer_Q, _transfer_T)

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


# ---- bankruptcy options (beyond the paper) -------------------------------
def test_partial_recovery_conserves_wealth_and_caps_realised_recovery():
    r = simulate(Params(model="L", delta_width=0.0, interest=0.3, recovery=0.5,
                        samples=(2 * 10**4,)), seed=2)
    assert r.defaults[0] > 0 and r.min_wealth[0] >= 0
    assert r.total_wealth[0] == pytest.approx(1000.0, rel=1e-9)  # haircut moves wealth, never destroys it
    assert r.recovery_rate[0] <= 0.5 + 1e-9


def test_plain_liquidation_recovers_less_than_the_full_claim():
    r = simulate(Params(model="L", samples=(10**5,)), seed=0)
    assert 0.5 < r.recovery_rate[0] < 1.0  # lender gets the borrower's remaining assets


def test_bankruptcy_cost_destroys_wealth():
    r = simulate(Params(model="L", delta_width=0.0, interest=0.3, bankruptcy_cost=0.3,
                        samples=(2 * 10**4,)), seed=2)
    assert r.total_wealth[0] < 1000.0 * 0.99


def test_partial_recovery_lowers_inequality_a_lot():
    base = g_at(10**5, model="L")
    assert g_at(10**5, model="L", recovery=0.75) < base - 0.15


def test_exclusion_reduces_defaults():
    free = simulate(Params(model="L", samples=(3 * 10**5,)), seed=1).defaults[0]
    banned = simulate(Params(model="L", exclusion_period=10**4, samples=(3 * 10**5,)), seed=1).defaults[0]
    assert banned < free


def test_bankruptcy_options_rejected_outside_bankrupt_L():
    with pytest.raises(ValueError):
        Params(model="J", recovery=0.5)
    with pytest.raises(ValueError):
        Params(model="L", insolvency="skip", bankruptcy_cost=0.1)
    with pytest.raises(ValueError):
        Params(model="L", recovery=1.5)


class _L:
    """Hand-built ledger state for unit tests of the helpers."""
    def __init__(self, n, debts, born_t=0, born=None):
        self.D = np.zeros((n, n)); self.born = np.zeros((n, n), dtype=np.int64)
        self.acc = np.zeros((n, n), dtype=np.int64)
        for (a, k), v in debts.items():
            self.D[a, k] = v; self.born[a, k] = (born or {}).get((a, k), born_t)
        self.dtot = self.D.sum(axis=1); self.st = np.zeros(3)
        self.qa = np.zeros(8 * n + 64, dtype=np.int64); self.qx = np.zeros(8 * n + 64)

    def garnish(self, m, before, g, t, term=0, rate=0.0, order=0, receipts=False):
        return _garnish_all(m, before, g, t, term, rate, order, receipts, self.D, self.born,
                            self.acc, self.dtot, self.st, self.qa, self.qx)


def test_garnish_all_pays_creditors_pro_rata_from_gains_only():
    before = np.ones(4)
    m = np.array([3.0, 0.5, 1.0, 2.0])           # gains: +2, loss, none, +1
    L = _L(4, {(0, 1): 3.0, (0, 2): 1.0, (1, 2): 10.0, (3, 0): 0.5})
    repaid, forgiven = L.garnish(m, before, 0.5, 1)
    D = L.D
    assert (repaid, forgiven) == (pytest.approx(1.0 + 0.5), 0.0)  # 0.5*2 and min(0.5*1, 0.5)
    assert m.sum() == pytest.approx(6.5)                          # pure transfer
    # agent 0 pays 1.0 split 3:1 to agents 1 and 2, agent 3 pays 0.5 to agent 0
    assert D[0, 1] == pytest.approx(3.0 - 0.75) and D[0, 2] == pytest.approx(1.0 - 0.25)
    assert D[1, 2] == 10.0 and D[3, 0] == pytest.approx(0.0)      # agent 1 had no gain
    assert m.tolist() == pytest.approx([3.0 - 1.0 + 0.5, 0.5 + 0.75, 1.0 + 0.25, 1.5])


def test_expired_debt_is_forgiven_not_paid():
    L = _L(3, {(0, 1): 2.0}, born_t=5)
    m = np.array([2.0, 1.0, 1.0])
    repaid, forgiven = L.garnish(m, np.ones(3), 1.0, 100, term=50)  # age 95 > 50
    assert repaid == 0.0 and forgiven == pytest.approx(2.0)
    assert L.D.sum() == 0.0 and m.tolist() == [2.0, 1.0, 1.0]


def test_seize_receivables_pro_rata_assignment_and_setoff():
    # defaulter 0 owes lender 1 an unpaid 3; 2 owes 0 six, 1 owes 0 two (set-off), 3 owes 0 four
    L = _L(4, {(2, 0): 6.0, (1, 0): 2.0, (3, 0): 4.0})
    seized, _ = _seize_receivables(0, 1, 3.0, 10, 0, 0.0, L.D, L.born, L.acc, L.dtot, L.st)
    D = L.D
    assert seized == pytest.approx(3.0)                  # f = 3/12 = 0.25 of every receivable
    assert D[2, 0] == pytest.approx(4.5) and D[2, 1] == pytest.approx(1.5)   # assigned to 1
    assert D[3, 0] == pytest.approx(3.0) and D[3, 1] == pytest.approx(1.0)
    assert D[1, 0] == pytest.approx(1.5) and D[1, 1] == 0.0                  # set off, not self-debt
    assert L.dtot[1] == pytest.approx(1.5) and L.dtot[2] == pytest.approx(6.0)   # 2's total unchanged


def test_repayment_order_fifo_lifo_pro_rata():
    debts = {(0, 1): 1.0, (0, 2): 1.0, (0, 3): 1.0}
    born = {(0, 1): 30, (0, 2): 10, (0, 3): 20}          # oldest default: creditor 2, then 3, then 1
    for order, expect in ((0, [1/3, 1/3, 1/3]), (1, [0.0, 1.0, 0.5]), (2, [1.0, 0.0, 0.5])):
        L = _L(4, debts, born=born)
        m = np.array([3.0, 1.0, 1.0, 1.0])               # debtor 0 gains 2, pays g*2 = 1.5
        L.garnish(m, np.ones(4), 0.75, 40, order=order)
        paid = (1.0 - L.D[0, 1:]).tolist()
        if order == 0:
            assert paid == pytest.approx([0.5, 0.5, 0.5])
        else:
            assert paid == pytest.approx([expect[0], expect[1], expect[2]])
        assert m.sum() == pytest.approx(6.0)


def test_interest_compounds_per_step_and_is_counted():
    L = _L(2, {(0, 1): 2.0}, born_t=0)
    m = np.array([2.0, 1.0])
    L.garnish(m, np.ones(2), 1.0, 100, rate=0.01)         # gain 1 -> pays min(1, grown debt)
    grown = 2.0 * 1.01 ** 100
    assert L.st[0] == pytest.approx(grown - 2.0)
    assert L.D[0, 1] == pytest.approx(grown - 1.0)        # one unit repaid out of the grown debt


def test_receipts_cascade_passes_on_to_the_recipients_creditors():
    # 0 owes 1 (3); 1 owes 2 (5). 0 gains 2 and pays g*2 = 2 to 1; with receipts, 1 passes
    # g*2 = 2 on to 2. Without receipts 1 keeps it.
    for receipts, final in ((False, [1.0, 3.0, 1.0]), (True, [1.0, 1.0, 3.0])):
        L = _L(3, {(0, 1): 3.0, (1, 2): 5.0})
        m = np.array([3.0, 1.0, 1.0])
        L.garnish(m, np.ones(3), 1.0, 10, receipts=receipts)
        assert m.tolist() == pytest.approx(final)
        assert L.st[1] == pytest.approx(2.0 if receipts else 0.0)


def test_ledger_conserves_wealth_and_repays_from_transfers():
    # dw=0 conserves wealth exactly; T transfers are the only income a defaulter has.
    r = simulate(Params(model="L", delta_width=0.0, interest=0.3, garnishment_rate=1.0,
                        redistribution="T", transfer_period=1000, samples=(2 * 10**4,)), seed=1)
    assert r.total_wealth[0] == pytest.approx(1000.0, rel=1e-9)
    assert r.debt_repaid[0] > 0 and r.debt_outstanding[0] > 0 and r.min_wealth[0] >= 0


def test_ledger_leaves_exchange_dynamics_unchanged_until_first_default_income():
    # No garnishment income exists before a defaulter's wealth rises, so g>0 must not
    # change the first 1e3 steps (no defaults yet).
    a = simulate(Params(model="L", samples=(1000,)), seed=5).wealth
    b = simulate(Params(model="L", garnishment_rate=0.5, samples=(1000,)), seed=5).wealth
    assert np.allclose(a, b)


def test_lender_exclusion_changes_default_count_and_never_deadlocks():
    base = simulate(Params(model="L", samples=(10**5,)), seed=1)
    lex = simulate(Params(model="L", lender_exclusion_period=10**4, samples=(10**5,)), seed=1)
    assert lex.defaults[0] != base.defaults[0]
    # With 2 agents, once one has defaulted and is banned from lending, every step in
    # which the other borrows has no eligible lender: bounded resampling -> idle step.
    tiny = simulate(Params(model="L", n_agents=2, interest=0.3, delta_width=0.3,
                           lender_exclusion_period=10**6, samples=(5000,)), seed=0)
    assert tiny.defaults[0] > 0 and tiny.idle_steps[0] > 0


def test_multi_debt_ledger_conserves_wealth_and_keeps_holdings_nonnegative():
    kw = dict(model="L", delta_width=0.0, interest=0.3, garnishment_rate=1.0,
              redistribution="T", transfer_period=1000, multi_debt=True,
              seize_receivables=True, samples=(2 * 10**4,))
    r = simulate(Params(**kw), seed=1)
    assert r.total_wealth[0] == pytest.approx(1000.0, rel=1e-9)
    assert r.min_wealth[0] >= 0 and r.debt_repaid[0] > 0 and r.claims_seized[0] > 0


def test_multi_debt_owes_more_than_single_debt():
    # Single mode drops a second default; multi mode records it, so open debt is larger.
    one = simulate(Params(model="L", garnishment_rate=0.5, samples=(10**5,)), seed=2)
    many = simulate(Params(model="L", garnishment_rate=0.5, multi_debt=True, samples=(10**5,)), seed=2)
    assert many.debt_outstanding[0] > one.debt_outstanding[0]


def test_debt_term_forgives_and_bounds_open_debt():
    base = simulate(Params(model="L", garnishment_rate=0.5, multi_debt=True, samples=(2 * 10**5,)), seed=2)
    short = simulate(Params(model="L", garnishment_rate=0.5, multi_debt=True, debt_term=2000,
                            samples=(2 * 10**5,)), seed=2)
    assert short.debt_forgiven[0] > 0 and base.debt_forgiven[0] == 0
    assert short.debt_outstanding[0] < 0.5 * base.debt_outstanding[0]


def test_material_defaults_are_a_subset_and_net_worth_share_is_a_share():
    r = simulate(Params(model="L", garnishment_rate=1.0, multi_debt=True, samples=(10**5,)), seed=1)
    assert 0 < r.material_defaults[0] < r.defaults[0]
    assert 0 < r.neg_net_worth_share[0] <= 1


def test_ledger_extensions_need_the_ledger():
    for kw in (dict(multi_debt=True), dict(debt_term=10), dict(seize_receivables=True)):
        with pytest.raises(ValueError):
            Params(model="L", **kw)
    with pytest.raises(ValueError):
        Params(model="L", garnishment_rate=0.5, seize_receivables=True)  # needs multi_debt


def test_ledger_rejects_combination_with_recovery():
    with pytest.raises(ValueError):
        Params(model="L", garnishment_rate=0.5, recovery=0.5)
    with pytest.raises(ValueError):
        Params(model="J", lender_exclusion_period=10)


def test_extended_ledger_options_conserve_wealth_end_to_end():
    kw = dict(model="L", delta_width=0.0, interest=0.3, garnishment_rate=1.0,
              redistribution="T", transfer_period=1000, multi_debt=True,
              seize_receivables=True, repayment_order="fifo", receipts_are_income=True,
              samples=(2 * 10**4,))
    r = simulate(Params(**kw), seed=1)
    assert r.total_wealth[0] == pytest.approx(1000.0, rel=1e-9)   # interest and priority never move holdings
    assert r.min_wealth[0] >= 0 and r.debt_repaid[0] > 0
    rate = simulate(Params(**dict(kw, debt_interest=1e-4)), seed=1)
    assert rate.total_wealth[0] == pytest.approx(1000.0, rel=1e-9)
    assert rate.interest_accrued[0] > 0 and rate.debt_outstanding[0] > r.debt_outstanding[0]


def test_cascade_never_leaves_negative_holdings():
    # Regression: roundoff in pass-through repayments used to leave holdings at -1e-15,
    # which a negative lender compounds by (1-lam)*rho per loan (min wealth reached -700
    # by t=3e5). A debtor can only repay what it holds.
    r = simulate(Params(model="L", garnishment_rate=1.0, multi_debt=True,
                        receipts_are_income=True, samples=(10**5, 3 * 10**5)), seed=1)
    assert (r.min_wealth >= 0).all()
    assert r.passed_on[-1] > 0 and r.total_wealth[-1] > 0


def test_new_ledger_options_validation():
    with pytest.raises(ValueError):
        Params(model="L", garnishment_rate=0.5, repayment_order="fifo")          # needs multi_debt
    with pytest.raises(ValueError):
        Params(model="L", debt_interest=1e-4)                                    # needs ledger
    with pytest.raises(ValueError):
        Params(model="L", garnishment_rate=0.5, repayment_order="newest")
    with pytest.raises(ValueError):
        Params(model="L", garnishment_rate=0.5, debt_interest=-1.0)


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
