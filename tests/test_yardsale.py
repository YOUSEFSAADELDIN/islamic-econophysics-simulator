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
