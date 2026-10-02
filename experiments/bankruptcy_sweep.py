"""Sweep the L-model bankruptcy options and score each against Kato's L figures.

    python experiments/bankruptcy_sweep.py [--seeds 10]

For every setting: Gini of L alone, of L-T (xi=0.5), defaults, total wealth and the
realised recovery rate, plus the mean/max abs deviation from the digitised Kato L
curves (Fig 4 rho=.05 and rho=0, Fig 5 tp=1e5, Fig 6 both tp, Fig 7).
"""

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from islamic_econophysics import kato2022_reference as ref  # noqa: E402
from islamic_econophysics.kato2022 import Params, simulate  # noqa: E402

N_SEEDS = 10
GRID = {
    "baseline (plain liquidation)": {},
    **{f"recovery={r}": dict(recovery=r) for r in (0.75, 0.5, 0.25, 0.0)},
    **{f"cost={c}": dict(bankruptcy_cost=c) for c in (0.1, 0.25, 0.5)},
    **{f"exclusion={e:g}": dict(exclusion_period=e) for e in (10**2, 10**3, 10**4)},
    **{f"lender exclusion={e:g}": dict(lender_exclusion_period=e) for e in (10**2, 10**3, 10**4)},
    **{f"ledger g={g}": dict(garnishment_rate=g) for g in (0.1, 0.5, 1.0)},
    # extended ledger (all at g=1; multi-creditor debts accumulate rather than being dropped)
    "ledger multi g=1": dict(garnishment_rate=1.0, multi_debt=True),
    "ledger multi+seize g=1": dict(garnishment_rate=1.0, multi_debt=True, seize_receivables=True),
    "ledger multi term=1e4 g=1": dict(garnishment_rate=1.0, multi_debt=True, debt_term=10**4),
    "ledger multi term=1e3 g=1": dict(garnishment_rate=1.0, multi_debt=True, debt_term=10**3),
    "ledger multi+seize+term=1e4 g=1": dict(garnishment_rate=1.0, multi_debt=True,
                                            seize_receivables=True, debt_term=10**4),
}


def run(args):
    kw, seed = args
    r = simulate(Params(**kw), seed=seed)
    return (r.gini, r.total_wealth, r.defaults, r.recovery_rate,
            r.debt_outstanding, r.debt_repaid, r.material_defaults, r.neg_net_worth_share,
            r.debt_forgiven, r.claims_seized)


def mean_over_seeds(pool, kw, seeds):
    out = pool.map(run, [(kw, s) for s in range(seeds)])
    return [np.nanmean([o[i] for o in out], axis=0) for i in range(10)]


def evaluate(pool, opts, seeds):
    base = dict(model="L", interest=0.05, delta_width=0.1, **opts)
    devs, row = [], {}
    # Fig 4, L rho=0.05 and rho=0 (exchange only)
    for rho, key in ((0.05, "L rho=0.05"), (0.0, "L rho=0")):
        g, tot, dfl, rec, dout, drep, mat, neg, forg, seiz = mean_over_seeds(
            pool, dict(base, interest=rho, samples=ref.FIG4_TIMES), seeds)
        devs.append(np.abs(g - ref.FIG4[key]))
        row[f"fig4 {key}"] = g.tolist()
        if rho == 0.05:
            row.update(defaults_1e6=float(dfl[-1]), total_per_agent_1e6=float(tot[-1] / 1000),
                       recovery_rate=float(rec[-1]) if np.isfinite(rec[-1]) else None,
                       debt_outstanding_per_agent_1e6=float(dout[-1] / 1000),
                       debt_repaid_1e6=float(drep[-1]),
                       material_defaults_1e6=float(mat[-1]), neg_net_share_1e6=float(neg[-1]),
                       debt_forgiven_1e6=float(forg[-1]), claims_seized_1e6=float(seiz[-1]))
    # Fig 5, L-T at tp=1e5
    g = mean_over_seeds(pool, dict(base, redistribution="T", transfer_rate=0.5,
                                   transfer_period=10**5, samples=ref.FIG4_TIMES), seeds)[0]
    devs.append(np.abs(g - ref.FIG5[("L", 10**5)])); row["fig5 L-T tp=1e5"] = g.tolist()
    # Fig 6 (xi sweep, t=1e6), both tp
    for tp in (10**4, 10**5):
        gs = [mean_over_seeds(pool, dict(base, redistribution="T", transfer_rate=xi,
                                         transfer_period=tp, samples=(10**6,)), seeds)[0][0]
              for xi in ref.FIG6_XI]
        devs.append(np.abs(np.array(gs) - ref.FIG6[("L", tp)])); row[f"fig6 tp={tp:g}"] = gs
    # Fig 7 (quantile sweep, t=1e6)
    gs = [mean_over_seeds(pool, dict(base, redistribution="Q", transfer_rate=0.5,
                                     transfer_period=10**5, quantile=float(q),
                                     samples=(10**6,)), seeds)[0][0] for q in ref.FIG7_Q]
    devs.append(np.abs(np.array(gs) - ref.FIG7["L"])); row["fig7"] = gs
    d = np.concatenate(devs)
    row.update(mean_dev=float(d.mean()), max_dev=float(d.max()),
               fig4_rho0_dev_1e6=float(devs[1][-1]))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=N_SEEDS)
    ap.add_argument("--out", default="results/bankruptcy_sweep")
    a = ap.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    res = {}
    with Pool() as pool:
        for name, opts in GRID.items():
            res[name] = evaluate(pool, opts, a.seeds)
            print(f"{name:32s} mean dev {res[name]['mean_dev']:.3f}  max {res[name]['max_dev']:.3f}", flush=True)
    Path(a.out, "sweep.json").write_text(json.dumps(res))
    lines = [f"# Bankruptcy options vs Kato's L figures ({a.seeds} seeds)\n",
             "Deviation = |sim - Kato| over the 7+7+7+7+7+10 = 45 digitised L points "
             "(Fig 4 rho=.05/0, Fig 5 tp=1e5, Fig 6 tp=1e4/1e5, Fig 7).\n",
             "| setting | mean dev | max dev | L alone t=1e6 (rho=.05) | L rho=0 t=1e6 | L-T tp=1e5 t=1e6 | defaults (all / material, 1e6) | wealth/agent | realised recovery | open debt/agent | repaid | forgiven | seized | agents with net worth < 0 |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n, r in res.items():
        rec = "-" if r["recovery_rate"] is None else f"{r['recovery_rate']:.2f}"
        lines.append(f"| {n} | {r['mean_dev']:.3f} | {r['max_dev']:.3f} | {r['fig4 L rho=0.05'][-1]:.2f} | "
                     f"{r['fig4 L rho=0'][-1]:.2f} | {r['fig5 L-T tp=1e5'][-1]:.2f} | "
                     f"{r['defaults_1e6']:.0f} / {r['material_defaults_1e6']:.0f} | {r['total_per_agent_1e6']:.2f} | {rec} | "
                     f"{r['debt_outstanding_per_agent_1e6']:.2f} | {r['debt_repaid_1e6']:.0f} | "
                     f"{r['debt_forgiven_1e6']:.0f} | {r['claims_seized_1e6']:.0f} | {r['neg_net_share_1e6']:.2f} |")
    Path(a.out, "sweep.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
