"""Zakat paid into a pool and split over the eight Qur'anic categories (Q 9:60).

    python experiments/zakat_categories.py [--parts A,B,C,D,E] [--seeds 6]

Base economy: N=2000, beta=0.1, waa=0.1, 20 trades per agent-year, 6000 lunar years,
zakat 2.5% on wealth above a nisab of 0.5 x mean, poverty line 0.2 x mean wealth.
  A  one category at a time (100% of the pool) and the stylised allocation profiles
  B  administrative overhead: weight of the amilin, rest to the poor and needy
  C  zakat rate vs poverty headcount for four profiles (rate needed to cut poverty)
  D  does splitting the pool change the critical rate (tail exponent = 1)?
  E  nisab sensitivity
Each part is cached in results/zakat_categories/part_X.json.
"""

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from islamic_econophysics.yardsale import CATEGORIES, PROFILES, YSParams, simulate  # noqa: E402
from zakat_sweep import hill_alpha, critical_rate  # noqa: E402

BASE = dict(n_agents=2000, beta=0.1, waa=0.1, trades_per_agent_year=20.0, years=6000, n_samples=40,
            zakat_rate=0.025, nisab=0.5, poverty_line=0.2, recipients="eight")
METRICS = ("gini", "poverty", "poverty_pre", "bottom50", "top10pct")


def one(args):
    kw, seed = args
    r = simulate(YSParams(**kw), seed=seed)
    d = {m: r.steady(m) for m in METRICS}
    d["alpha"] = hill_alpha(r.wealth)
    st = r.stats()
    d["to_poorest20"] = st.get("share_to_poorest_20pct", float("nan"))
    d["unspent"] = st.get("share_unspent_returned", float("nan"))
    for c in CATEGORIES:
        d["share_" + c] = st.get("share_" + c, float("nan"))
    return d


def cfg(pool, kw, seeds):
    out = pool.map(one, [(dict(BASE, **kw), s) for s in range(seeds)])
    return {k: [float(np.nanmean([o[k] for o in out])), float(np.nanstd([o[k] for o in out]))] for k in out[0]}


def single(k):
    w = [0.0] * 8
    w[k] = 1.0
    return tuple(w)


def part_a(pool, seeds):
    res = {"no zakat": cfg(pool, dict(zakat_rate=0.0), seeds)}
    for name, wts in PROFILES.items():
        res[f"profile:{name}"] = cfg(pool, dict(category_weights=wts), seeds)
    for k, c in enumerate(CATEGORIES):
        res[f"only:{c}"] = cfg(pool, dict(category_weights=single(k)), seeds)
        print("  A", c, round(res[f"only:{c}"]["poverty_pre"][0], 3), flush=True)
    return res


def part_b(pool, seeds):
    res = {}
    for a in (0.0, 0.125, 0.25, 0.4):
        wts = (0.5 * (1 - a), 0.5 * (1 - a), a, 0, 0, 0, 0, 0)
        res[f"amilin={a}"] = cfg(pool, dict(category_weights=wts), seeds)
    return res


def part_c(pool, seeds):
    res = {}
    for name in ("all_public", "poor_only", "needy_first", "eight_equal"):
        res[name] = {str(z): cfg(pool, dict(category_weights=PROFILES[name], zakat_rate=z), seeds)
                     for z in (0.005, 0.01, 0.025, 0.05, 0.1)}
    return res


def part_d(pool, seeds):
    zs = [0.0005, 0.001, 0.0015, 0.0025, 0.004, 0.006, 0.01, 0.02]
    res = {}
    for name in ("all_public", "poor_only", "eight_equal"):
        rows = [cfg(pool, dict(category_weights=PROFILES[name], zakat_rate=z), seeds) for z in zs]
        res[name] = dict(z=zs, alpha=[r["alpha"][0] for r in rows],
                         z_c=critical_rate(zs, [r["alpha"][0] for r in rows]))
        print("  D", name, res[name]["z_c"], flush=True)
    return res


def part_e(pool, seeds):
    res = {}
    for nis in (0.3, 0.5, 1.0):
        for name in ("poor_only", "eight_equal"):
            res[f"nisab={nis},{name}"] = cfg(pool, dict(category_weights=PROFILES[name], nisab=nis), seeds)
    return res


def write_tables(out, R):
    L = []
    if "A" in R:
        A = R["A"]
        L += ["## A. One category at a time, and stylised profiles (z=2.5%, nisab 0.5, N=2000, waa=0.1)\n",
              "| allocation | Gini | poverty before zakat / after (<0.2 mean) | bottom-50% share | top-10% share | pool reaching poorest 20% | unspent, returned |",
              "|---|---|---|---|---|---|---|"]
        for k, v in A.items():
            f = lambda m, d=3: f"{v[m][0]:.{d}f}"
            to = "-" if np.isnan(v["to_poorest20"][0]) else f"{v['to_poorest20'][0]:.2f}"
            un = "-" if np.isnan(v["unspent"][0]) else f"{v['unspent'][0]:.2f}"
            L.append(f"| {k} | {f('gini', 2)} | {f('poverty_pre')} / {f('poverty')} | {f('bottom50')} | {f('top10pct', 2)} | {to} | {un} |")
        L.append("")
    if "B" in R:
        L += ["## B. Administrative overhead (weight of the amilin; the rest split between poor and needy)\n",
              "| amilin weight | Gini | poverty before zakat | bottom-50% share | pool reaching poorest 20% |", "|---|---|---|---|---|"]
        for k, v in R["B"].items():
            L.append(f"| {k.split('=')[1]} | {v['gini'][0]:.2f} | {v['poverty_pre'][0]:.3f} | {v['bottom50'][0]:.3f} | {v['to_poorest20'][0]:.2f} |")
        L.append("")
    if "C" in R:
        L += ["## C. Zakat rate vs poverty headcount before zakat arrives (Gini in brackets)\n", "| zakat/yr | " + " | ".join(R["C"].keys()) + " |",
              "|---|" + "---|" * len(R["C"])]
        for z in next(iter(R["C"].values())):
            L.append(f"| {float(z):.3f} | " + " | ".join(f"{R['C'][n][z]['poverty_pre'][0]:.3f} (Gini {R['C'][n][z]['gini'][0]:.2f})" for n in R["C"]) + " |")
        L.append("")
    if "D" in R:
        L += ["## D. Critical rate (tail exponent = 1) by allocation\n", "| allocation | z_c |", "|---|---|"]
        for n, v in R["D"].items():
            L.append(f"| {n} | {v['z_c']:.5f} |" if v["z_c"] else f"| {n} | n/a |")
        L.append("")
    if "E" in R:
        L += ["## E. Nisab sensitivity (poverty line fixed at 0.2 x mean)\n",
              "| setting | Gini | poverty before zakat | bottom-50% share |", "|---|---|---|---|"]
        for k, v in R["E"].items():
            L.append(f"| {k} | {v['gini'][0]:.2f} | {v['poverty_pre'][0]:.3f} | {v['bottom50'][0]:.3f} |")
    (out / "tables.md").write_text("\n".join(L) + "\n")


def figures(out, R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if "A" in R:
        A = R["A"]
        names = [k for k in A if k != "no zakat"]
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        for a, key, ttl in zip(ax, ("poverty_pre", "gini"), ("poverty headcount before zakat (< 0.2 x mean)", "Gini")):
            vals = [A[k][key][0] for k in names]
            errs = [A[k][key][1] for k in names]
            a.barh(names, vals, xerr=errs, color=["C0" if k.startswith("profile") else "C1" for k in names])
            a.axvline(A["no zakat"][key][0], color="k", ls="--", lw=.8)
            a.set_title(ttl + " (dashed = no zakat)")
        fig.tight_layout(); fig.savefig(out / "league_table.png", dpi=130); plt.close(fig)
    if "C" in R:
        fig, ax = plt.subplots(figsize=(6, 4.2))
        for n, v in R["C"].items():
            zs = [float(z) for z in v]
            ax.plot(zs, [v[z]["poverty_pre"][0] for z in v], "o-", label=n)
        ax.set_xscale("log"); ax.set_xlabel("zakat rate per lunar year"); ax.set_ylabel("poverty headcount before zakat")
        ax.legend(); fig.tight_layout(); fig.savefig(out / "poverty_vs_rate.png", dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="A,B,C,D,E")
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--out", default="results/zakat_categories")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    R, fn = {}, dict(A=part_a, B=part_b, C=part_c, D=part_d, E=part_e)
    with Pool() as pool:
        for part in a.parts.split(","):
            f = out / f"part_{part}.json"
            if f.exists():
                R[part] = json.loads(f.read_text())
                continue
            print("part", part, flush=True)
            R[part] = fn[part](pool, a.seeds)
            f.write_text(json.dumps(R[part]))
    write_tables(out, R)
    figures(out, R)
    print("done")


if __name__ == "__main__":
    main()
