"""Reproduce Figs 3-7 of Kato (2022) and compare with the digitised figures.

    python experiments/replicate_kato2022.py [--seeds 20] [--outdir results/kato2022]

Writes results.json, comparison.md and fig3..fig7 PNGs. Each Gini value is the
mean over `--seeds` independent runs (the paper shows a single run, which is
why the error bars here are worth looking at).
"""

from __future__ import annotations

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from islamic_econophysics import kato2022_reference as ref  # noqa: E402
from islamic_econophysics.kato2022 import Params, simulate  # noqa: E402

T_FINAL = 10**6
# label -> model kwargs. "J0.2" etc. match the keys in kato2022_reference.
EXCHANGE = {
    "R": dict(model="R"),
    "L rho=0.05": dict(model="L", interest=0.05, delta_width=0.1),
    "L rho=0": dict(model="L", interest=0.0, delta_width=0.1),
    "J dw=0.2": dict(model="J", delta_width=0.2),
    "J dw=0.1": dict(model="J", delta_width=0.1),
}
SHORT = {"R": "R", "L": "L rho=0.05", "J0.2": "J dw=0.2", "J0.1": "J dw=0.1"}


def _one(args):
    kw, seed = args
    r = simulate(Params(**kw), seed=seed)
    return r.gini, r.total_wealth


def ensemble(pool, kw, seeds):
    out = pool.map(_one, [(kw, s) for s in range(seeds)])
    g = np.array([o[0] for o in out])
    ensemble.last_total = np.mean([o[1] for o in out], axis=0)  # diagnostic
    return g.mean(0), g.std(0), out


def compare(name, sim, target, tol=ref.TOL):
    dev = np.abs(np.asarray(sim) - np.asarray(target))
    return dict(series=name, max_dev=float(dev.max()),
                mean_dev=float(dev.mean()), ok=bool(dev.max() <= tol))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--outdir", default="results/kato2022")
    ap.add_argument("--insolvency", default="bankrupt",
                    help="L-model rule for a borrower who would go below 0 (see Params)")
    a = ap.parse_args()
    for kw in EXCHANGE.values():
        if kw["model"] == "L":
            kw["insolvency"] = a.insolvency
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    res, rows = {}, []

    with Pool() as pool:
        # ---- Fig 4: exchange only --------------------------------------
        fig4 = {}
        for name, kw in EXCHANGE.items():
            m, s, _ = ensemble(pool, dict(kw, samples=ref.FIG4_TIMES), a.seeds)
            fig4[name] = dict(mean=m.tolist(), sd=s.tolist(),
                              total_wealth=ensemble.last_total.tolist())
            rows.append(compare(f"Fig4 {name}", m, ref.FIG4[name]))
        res["fig4"] = fig4

        # ---- Fig 3: wealth histograms (single run, like the paper) ----
        fig3 = {}
        snap_t = (10**3, 10**4, 10**6)
        for name in ("R", "L rho=0.05", "J dw=0.1", "J dw=0.2"):
            r = simulate(Params(**EXCHANGE[name], samples=snap_t), seed=0)
            fig3[name] = {str(t): r.wealth[k].tolist() for k, t in enumerate(snap_t)}
        res["fig3"] = fig3

        # ---- Fig 5: T model over time ---------------------------------
        fig5 = {}
        tkw = dict(redistribution="T", transfer_rate=0.5)
        for tp in (10**4, 10**5):
            for key, name in SHORT.items():
                kw = dict(EXCHANGE[name], **tkw, transfer_period=tp,
                          samples=ref.FIG4_TIMES)
                m, s, _ = ensemble(pool, kw, a.seeds)
                fig5[f"{key}|{tp}"] = dict(mean=m.tolist(), sd=s.tolist())
                if (key, tp) in ref.FIG5:
                    rows.append(compare(f"Fig5 {key} tp={tp:g}", m, ref.FIG5[(key, tp)]))
                elif tp == 10**4:
                    plateau = float(np.mean(m[-3:]))
                    rows.append(dict(series=f"Fig5 {key} tp=1e4 plateau",
                                     max_dev=abs(plateau - ref.FIG5_PLATEAU_TP1E4[key]),
                                     mean_dev=abs(plateau - ref.FIG5_PLATEAU_TP1E4[key]),
                                     ok=abs(plateau - ref.FIG5_PLATEAU_TP1E4[key]) <= ref.TOL))
        res["fig5"] = fig5

        # ---- Fig 6: Gini vs transfer rate ------------------------------
        fig6 = {}
        for tp in (10**4, 10**5):
            for key, name in SHORT.items():
                ms, ss = [], []
                for xi in ref.FIG6_XI:
                    kw = dict(EXCHANGE[name], redistribution="T", transfer_rate=xi,
                              transfer_period=tp, samples=(T_FINAL,))
                    m, s, _ = ensemble(pool, kw, a.seeds)
                    ms.append(float(m[0])); ss.append(float(s[0]))
                fig6[f"{key}|{tp}"] = dict(mean=ms, sd=ss)
                if (key, tp) in ref.FIG6:
                    rows.append(compare(f"Fig6 {key} tp={tp:g}", ms, ref.FIG6[(key, tp)]))
        res["fig6"] = fig6

        # ---- Fig 7: quantile model --------------------------------------
        fig7 = {}
        for key in ("L", "J0.2", "J0.1"):
            ms, ss = [], []
            for q in ref.FIG7_Q:
                kw = dict(SHORT_KW(key), redistribution="Q", transfer_rate=0.5,
                          transfer_period=10**5, quantile=float(q), samples=(T_FINAL,))
                m, s, _ = ensemble(pool, kw, a.seeds)
                ms.append(float(m[0])); ss.append(float(s[0]))
            fig7[key] = dict(mean=ms, sd=ss)
            rows.append(compare(f"Fig7 {key}", ms, ref.FIG7[key]))
        res["fig7"] = fig7

    res["meta"] = dict(seeds=a.seeds, t_final=T_FINAL, tolerance=ref.TOL,
                       insolvency=a.insolvency)
    (out / "results.json").write_text(json.dumps(res))
    write_comparison(out / "comparison.md", rows, a.seeds, a.insolvency)
    make_figures(res, out)
    bad = [r["series"] for r in rows if not r["ok"]]
    print(f"{len(rows) - len(bad)}/{len(rows)} series within +/-{ref.TOL}")
    for b in bad:
        print("  outside tolerance:", b)


def SHORT_KW(key):
    return EXCHANGE[SHORT[key]]


def write_comparison(path, rows, seeds, insolvency):
    lines = [f"# Replication vs. digitised Kato (2022) figures ({seeds} seeds/point)\n",
             f"L-model insolvency rule: `{insolvency}`.\n",
             f"Tolerance +/-{ref.TOL} on the Gini index (figure reading error + Monte-Carlo noise).\n",
             "| series | max abs dev | mean abs dev | within tol |", "|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['series']} | {r['max_dev']:.3f} | {r['mean_dev']:.3f} | "
                     f"{'yes' if r['ok'] else '**NO**'} |")
    path.write_text("\n".join(lines) + "\n")


def make_figures(res, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"R": "#1b7837", "L rho=0.05": "#762a83", "L rho=0": "#c2a5cf",
              "J dw=0.2": "#2166ac", "J dw=0.1": "#67a9cf"}
    t = np.array(ref.FIG4_TIMES)

    # Fig 3
    fig, axs = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
    for ax, (name, snaps) in zip(axs.ravel(), res["fig3"].items()):
        for tt, w in snaps.items():
            w = np.asarray(w)
            ax.hist(w / w.mean(), bins=np.linspace(0, 3, 21), alpha=.55,
                    label=f"t=10^{int(np.log10(int(tt)))}")
        ax.set_title(name); ax.set_xlabel("wealth m / mean wealth"); ax.set_ylabel("frequency"); ax.legend()
    fig.suptitle("Fig 3 replication: wealth distributions (N=1000, one run)")
    fig.tight_layout(); fig.savefig(out / "fig3.png", dpi=130); plt.close(fig)

    # Fig 4
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for name, d in res["fig4"].items():
        ax.errorbar(t, d["mean"], d["sd"], color=colors[name], marker="o", ms=4, label=f"{name} (sim)")
        ax.plot(t, ref.FIG4[name], color=colors[name], ls=":", marker="x", ms=5, alpha=.8)
    ax.set_xscale("log"); ax.set_xlabel("time t"); ax.set_ylabel("Gini")
    ax.set_title("Fig 4: exchange only (dotted x = Kato, read off plot)")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(out / "fig4.png", dpi=130); plt.close(fig)

    # Fig 5
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, tp in zip(axs, (10**4, 10**5)):
        for key, name in SHORT.items():
            d = res["fig5"][f"{key}|{tp}"]
            ax.errorbar(t, d["mean"], d["sd"], color=colors[name], marker="o", ms=4, label=name)
            if (key, tp) in ref.FIG5 and key != "R":
                ax.plot(t, ref.FIG5[(key, tp)], color=colors[name], ls=":", marker="x")
        ax.set_xscale("log"); ax.set_title(f"T model, xi=0.5, tp=1e{int(np.log10(tp))}")
        ax.set_xlabel("time t")
    axs[0].set_ylabel("Gini"); axs[0].legend(fontsize=7)
    fig.suptitle("Fig 5 replication (dotted x = Kato)"); fig.tight_layout()
    fig.savefig(out / "fig5.png", dpi=130); plt.close(fig)

    # Fig 6
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for tp, ls in ((10**5, "-"), (10**4, "--")):
        for key, name in SHORT.items():
            d = res["fig6"][f"{key}|{tp}"]
            ax.errorbar(ref.FIG6_XI, d["mean"], d["sd"], color=colors[name], ls=ls, marker="o", ms=4,
                        label=f"{name}, tp=1e{int(np.log10(tp))}")
            if (key, tp) in ref.FIG6:
                ax.plot(ref.FIG6_XI, ref.FIG6[(key, tp)], color=colors[name], ls=":", marker="x", alpha=.6)
    ax.set_xlabel("transfer rate xi"); ax.set_ylabel("Gini (t=1e6)")
    ax.set_title("Fig 6 replication (dotted x = Kato)"); ax.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(out / "fig6.png", dpi=130); plt.close(fig)

    # Fig 7
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for key in ("L", "J0.2", "J0.1"):
        name = SHORT[key]; d = res["fig7"][key]
        ax.errorbar(ref.FIG7_Q, d["mean"], d["sd"], color=colors[name], marker="o", ms=4, label=f"{name}-Q")
        ax.plot(ref.FIG7_Q, ref.FIG7[key], color=colors[name], ls=":", marker="x", alpha=.7)
        ax.axhline(ref.FIG7_T_LEVELS[key], color=colors[name], ls="--", lw=.8)
    ax.set_xlabel("q (threshold = m_max / q)"); ax.set_ylabel("Gini (t=1e6)")
    ax.set_title("Fig 7 replication (dashed = T model, dotted x = Kato)")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(out / "fig7.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
