import numpy as np
import pytest

from islamic_econophysics.yardsale import YSParams, _zakat, simulate


def test_total_wealth_is_conserved_and_positive():
    for kw in (dict(), dict(waa=0.3), dict(zakat_rate=0.05, nisab=0.5, recipients="poor"),
               dict(waa=0.2, zakat_rate=0.025)):
        r = simulate(YSParams(n_agents=200, years=300, n_samples=5, **kw), seed=1)
        assert r.wealth.sum() == pytest.approx(200.0, rel=1e-9)
        assert (r.wealth > 0).all()


def test_zakat_levy_rules():
    w = np.array([0.2, 0.9, 1.0, 3.9])           # mean 1.5; nisab 0.5 -> threshold 0.75
    _zakat(w, 0.1, 0.5, 0)                        # all recipients
    assert w.sum() == pytest.approx(6.0)
    # 0.2 is below the nisab and pays nothing; the others pay 10% and the pool is shared equally
    pool = 0.1 * (0.9 + 1.0 + 3.9)
    assert w.tolist() == pytest.approx([0.2 + pool / 4, 0.81 + pool / 4, 0.9 + pool / 4, 3.51 + pool / 4])
    w = np.array([0.2, 0.9, 1.0, 3.9])
    _zakat(w, 0.1, 0.5, 1)                        # only the agent below the nisab receives
    assert w.sum() == pytest.approx(6.0)
    assert w.tolist() == pytest.approx([0.2 + pool, 0.81, 0.9, 3.51])


def test_nisab_above_everyone_means_no_zakat():
    a = simulate(YSParams(n_agents=100, years=200, n_samples=4, waa=0.1), seed=3)
    b = simulate(YSParams(n_agents=100, years=200, n_samples=4, waa=0.1, zakat_rate=0.2, nisab=1e6), seed=3)
    assert np.allclose(a.wealth, b.wealth)


def test_fair_game_mean_wealth_of_an_agent_is_a_martingale():
    # waa=0: E[w_i(t)] = 1 for every agent; average over many runs of a short game
    ws = np.mean([simulate(YSParams(n_agents=20, years=5, n_samples=1), seed=s).wealth for s in range(400)], axis=0)
    assert np.abs(ws - 1.0).max() < 0.1


def test_waa_favours_the_richer_agent():
    # same stream of coins' stakes; with waa the early leaders pull ahead: higher Gini after a short game
    g0 = np.mean([simulate(YSParams(n_agents=200, years=40, n_samples=1, waa=0.0), seed=s).gini[-1] for s in range(6)])
    g1 = np.mean([simulate(YSParams(n_agents=200, years=40, n_samples=1, waa=0.5), seed=s).gini[-1] for s in range(6)])
    assert g1 > g0 + 0.02


def test_zakat_lowers_inequality_monotonically():
    gs = [np.mean([simulate(YSParams(n_agents=300, years=1500, n_samples=10, waa=0.1, zakat_rate=z), seed=s).steady("gini")
                   for s in range(3)]) for z in (0.0, 0.005, 0.025, 0.1)]
    assert all(a > b for a, b in zip(gs, gs[1:]))


def test_validation():
    with pytest.raises(ValueError):
        YSParams(waa=1.5)
    with pytest.raises(ValueError):
        YSParams(recipients="rich")
    with pytest.raises(ValueError):
        YSParams(beta=1.0)


# ---- eight recipient categories ---------------------------------------------------
from islamic_econophysics.yardsale import CATEGORIES, PROFILES, _zakat_eight  # noqa: E402


def _state(w):
    w = np.array(w, dtype=float)
    return w, w.copy(), np.zeros(len(w), dtype=np.int64), np.zeros(11)


CAT = np.array([0.2, 3, 0.5, 0.2, 0.2, 0.2])  # poverty_line, riqab_years, gharimin_loss, fractions


def test_eight_conserves_wealth_for_every_profile_and_pays_out_the_whole_pool():
    for name, wts in PROFILES.items():
        rng = np.random.default_rng(1)
        w, wp, pov, rcv = _state(rng.pareto(1.5, size=50) + 0.01)
        w *= 50 / w.sum(); wp = w.copy()
        total = w.sum()
        pool = _zakat_eight(w, wp, pov, 0.1, 0.5, np.array(wts, float), CAT, rcv)
        assert w.sum() == pytest.approx(total), name
        assert (w > 0).all()
        assert rcv[8] == pytest.approx(pool)
        assert rcv[:8].sum() + rcv[10] == pytest.approx(pool)       # paid + returned = collected


def test_poor_and_needy_are_paid_by_gap_and_never_beyond_the_nisab():
    # mean = 1 (sum 6 over 6 agents); nisab 0.5 -> thr 0.5; poverty line 0.2
    w, wp, pov, rcv = _state([0.05, 0.1, 0.3, 0.4, 3.15, 2.0])
    pool = _zakat_eight(w, wp, pov, 0.1, 0.5, np.array(PROFILES["poor_only"], float), CAT, rcv)
    # payers: the two agents at/above 0.5 (3.15 and 2.0) -> pool = 0.515
    assert pool == pytest.approx(0.515)
    # fuqara get 0.2575 split by gap to 0.5 among agents < 0.2 (gaps .45, .4); masakin the rest among .3, .4
    assert w[0] == pytest.approx(0.05 + 0.2575 * 0.45 / 0.85) and w[1] == pytest.approx(0.1 + 0.2575 * 0.4 / 0.85)
    assert w[2] == pytest.approx(0.3 + 0.2575 * 0.2 / 0.3) and w[3] == pytest.approx(0.4 + 0.2575 * 0.1 / 0.3)
    assert w[2] <= 0.5 and w[3] <= 0.5


def test_unspent_category_share_is_returned_to_everyone():
    # nobody is stranded below the poverty line, so the fuqara share cannot be spent
    w, wp, pov, rcv = _state([1.0, 1.0, 1.0, 1.0])
    _zakat_eight(w, wp, pov, 0.1, 0.5, np.array([1, 0, 0, 0, 0, 0, 0, 0], float), CAT, rcv)
    assert rcv[0] == 0.0 and rcv[10] == pytest.approx(rcv[8])
    assert w.tolist() == pytest.approx([1.0] * 4)  # levy returned equally


def test_riqab_pays_only_persistent_poor_and_gharimin_only_big_losers():
    w, wp, pov, rcv = _state([0.05, 0.05, 0.4, 3.5])
    pov[:] = [5, 1, 0, 0]                                   # only agent 0 has been poor >= 3 years
    _zakat_eight(w, wp, pov, 0.1, 0.5, np.array([0, 0, 0, 0, 1, 0, 0, 0], float), CAT, rcv)
    assert rcv[4] > 0 and w[0] > 0.05 and w[1] == pytest.approx(0.05)
    w, wp, pov, rcv = _state([0.3, 0.45, 0.4, 2.85])
    wp = np.array([1.0, 0.5, 0.4, 2.9])                     # agent 0 lost 70% (> 50%), agent 1 only 10%
    _zakat_eight(w, wp, pov, 0.1, 0.5, np.array([0, 0, 0, 0, 0, 1, 0, 0], float), CAT, rcv)
    assert w[0] > 0.3 and w[1] == pytest.approx(0.45)
    assert w[0] - 0.3 <= 0.7 + 1e-12                        # never more than the loss


def test_eight_simulation_runs_conserves_and_is_reproducible():
    kw = dict(n_agents=200, years=200, n_samples=5, waa=0.1, zakat_rate=0.05, nisab=0.5, recipients="eight")
    a = simulate(YSParams(**kw), seed=2); b = simulate(YSParams(**kw), seed=2)
    assert np.array_equal(a.wealth, b.wealth)
    assert a.wealth.sum() == pytest.approx(200.0)
    st = a.stats()
    assert sum(st[f"share_{c}"] for c in CATEGORIES) + st["share_unspent_returned"] == pytest.approx(1.0)


def test_eight_validation():
    with pytest.raises(ValueError):
        YSParams(recipients="eight", nisab=0.0)
    with pytest.raises(ValueError):
        YSParams(recipients="eight", nisab=0.5, poverty_line=0.6)
    with pytest.raises(ValueError):
        YSParams(recipients="eight", nisab=0.5, category_weights=(0,) * 8)


def test_poverty_before_zakat_is_at_least_poverty_after_when_the_poor_are_paid():
    r = simulate(YSParams(n_agents=300, years=1500, n_samples=20, waa=0.1, zakat_rate=0.05, nisab=0.5,
                          recipients="eight", category_weights=PROFILES["poor_only"]), seed=1)
    assert r.steady("poverty_pre") >= r.steady("poverty") - 1e-12
    assert r.steady("poverty_pre") > 0.0
