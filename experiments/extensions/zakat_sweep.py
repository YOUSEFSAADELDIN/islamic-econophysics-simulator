"""Zakat-rate sweep on the yard-sale model with wealth-attained advantage (WAA).

    python experiments/extensions/zakat_sweep.py [--parts A,B,C,D] [--seeds 6] [--out results/extensions/zakat]

Parts (each cached in results/extensions/zakat/part_X.json; delete to recompute):
  A  dose-response: waa x zakat rate -> Gini, top shares, Pareto tail exponent, critical rate
  B  how the critical rate scales with trades per agent-year and with the stake beta
  C  nisab and recipients (all vs poor) at a fixed waa and the 2.5% rate
  D  population size: is the critical rate a property of the model or of N?

Critical rate: the zakat rate at which the Pareto exponent of the stationary wealth tail
(Hill estimator on the top 5% of agents) crosses 1, i.e. below which the mean is carried by
the extreme tail and a few agents hold a finite share however large N is. Interpolated
log-linearly in z between the grid points that bracket 1.
"""

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from islamic_econophysics.extensions.yardsale import YSParams, simulate  # noqa: E402

YEARS = 6000
Z_GRID = [0.0005, 0.001, 0.0015, 0.0025, 0.004, 0.006, 0.01, 0.015, 0.025, 0.04, 0.07, 0.1]
WAA_GRID = [0.0, 0.02, 0.05, 0.1, 0.2]


def hill_alpha(w, frac=0.05):
    w = np.sort(w)[::-1]
    k = max(20, int(len(w) * frac))
    return float(1.0 / np.mean(np.log(w[:k] / w[k])))


def one(args):
    kw, seed = args
    p = YSParams(years=kw.pop("years", YEARS), n_samples=40, **kw)
    r = simulate(p, seed=seed)
    w = np.sort(r.wealth)[::-1]
    n = len(w)
    return dict(gini=r.steady("gini"), top10=r.steady("top10pct"), top1=r.steady("top1"),
                top1pct=float(w[:max(1, n // 100)].sum() / w.sum()),
                bottom50=float(w[n // 2:].sum() / w.sum()), alpha=hill_alpha(w))


def cfg(pool, kw, seeds):
    out = pool.map(one, [(dict(kw), s) for s in range(seeds)])
    return {k: [float(np.mean([o[k] for o in out])), float(np.std([o[k] for o in out]))] for k in out[0]}


def critical_rate(zs, alphas):
    """z where alpha crosses 1 (log-linear interpolation); None if it never does."""
    a = np.asarray(alphas)
    for i in range(len(zs) - 1):
        if a[i] <= 1.0 < a[i + 1]:
            t = (1.0 - a[i]) / (a[i + 1] - a[i])
            return float(np.exp(np.log(zs[i]) + t * (np.log(zs[i + 1]) - np.log(zs[i]))))
    return None if a[0] > 1.0 or a[-1] <= 1.0 else None


def curve(pool, base, zs, seeds):
    rows = [cfg(pool, dict(base, zakat_rate=z), seeds) for z in zs]
    return dict(z=list(zs), rows=rows, z_c=critical_rate(zs, [r["alpha"][0] for r in rows]))


def part_a(pool, seeds):
    return {str(w): curve(pool, dict(n_agents=2000, waa=w), Z_GRID, seeds) for w in WAA_GRID}


def part_b(pool, seeds):
    out = {}
    for ny in (5, 20, 80):
        for beta in (0.05, 0.1, 0.2):
            if ny == 80 and beta != 0.1:
                continue
            for w in (0.05, 0.2):
                zs = [z for z in Z_GRID if z <= 0.1]
                out[f"ny={ny},beta={beta},waa={w}"] = curve(
                    pool, dict(n_agents=1000, waa=w, beta=beta, trades_per_agent_year=ny), zs, seeds)
                print("  B", f"ny={ny} beta={beta} waa={w}", out[f"ny={ny},beta={beta},waa={w}"]["z_c"], flush=True)
    return out


def part_c(pool, seeds):
    out = {}
    for rec in ("all", "poor"):
        for nis in (0.0, 0.1, 0.5, 1.0):
            out[f"{rec},nisab={nis}"] = {
                f"waa={w}": cfg(pool, dict(n_agents=2000, waa=w, zakat_rate=0.025, nisab=nis, recipients=rec), seeds)
                for w in (0.05, 0.2)}
            # critical rate with this nisab/recipient rule
            out[f"{rec},nisab={nis}"]["z_c(waa=0.1)"] = curve(
                pool, dict(n_agents=2000, waa=0.1, nisab=nis, recipients=rec), Z_GRID, seeds)["z_c"]
            print("  C", rec, nis, out[f"{rec},nisab={nis}"]["z_c(waa=0.1)"], flush=True)
    return out


def part_d(pool, seeds):
    out = {}
    for n in (500, 2000, 8000):
        out[str(n)] = curve(pool, dict(n_agents=n, waa=0.1), Z_GRID, seeds)
        print("  D", n, out[str(n)]["z_c"], flush=True)
    return out


def make_figures(out, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    A = res.get("A")
    if A:
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        cols = plt.cm.viridis(np.linspace(0, 0.9, len(A)))
        for c, (w, d) in zip(cols, A.items()):
            z = np.array(d["z"])
            for a, key in zip(ax, ("gini", "top10", "alpha")):
                a.plot(z, [r[key][0] for r in d["rows"]], "o-", ms=3, color=c, label=f"waa={w}")
        ax[0].set_ylabel("Gini"); ax[1].set_ylabel("top-10% wealth share"); ax[2].set_ylabel("Pareto exponent of the tail")
        ax[2].axhline(1, color="k", ls="--", lw=.8)
        for a in ax:
            a.set_xscale("log"); a.set_xlabel("zakat rate per lunar year"); a.axvline(0.025, color="r", lw=.8, ls=":")
        ax[0].legend(fontsize=7)
        fig.suptitle("Yard sale + WAA + zakat (N=2000): dose-response; dotted red = 2.5%; alpha=1 is the oligarchy threshold")
        fig.tight_layout(); fig.savefig(out / "dose_response.png", dpi=130); plt.close(fig)
    B = res.get("B")
    if B:
        coef, rows = fit_critical(res)
        fig, ax = plt.subplots(figsize=(5.2, 5))
        pred = np.array([coef[0] * r[0] * r[1] * r[2] + coef[1] * r[0] * r[1] ** 2 for r in rows])
        meas = np.array([r[3] for r in rows])
        ax.loglog(pred, meas, "o")
        lim = [min(pred.min(), meas.min()) * 0.8, max(pred.max(), meas.max()) * 1.2]
        ax.plot(lim, lim, "k--", lw=.8)
        ax.axhline(0.025, color="r", ls=":", lw=.8)
        ax.set_xlabel(f"fit  {coef[0]:.4f} n_y beta waa + {coef[1]:.4f} n_y beta^2")
        ax.set_ylabel("measured critical zakat rate z_c (per lunar year)")
        ax.set_title("critical rate collapses onto the fit; red dotted = 2.5%")
        fig.tight_layout(); fig.savefig(out / "critical_scaling.png", dpi=130); plt.close(fig)


def fit_critical(res):
    """Least-squares z_c = a n_y beta waa + b n_y beta^2 over all (n_y, beta, waa) points of parts A and B."""
    rows = []
    for k, d in res["B"].items():
        if d["z_c"] is not None:
            rows.append((float(k.split("ny=")[1].split(",")[0]), float(k.split("beta=")[1].split(",")[0]),
                         float(k.split("waa=")[1]), d["z_c"]))
    for w, d in res["A"].items():
        if d["z_c"] is not None:
            rows.append((20.0, 0.1, float(w), d["z_c"]))
    X = np.array([[r[0] * r[1] * r[2], r[0] * r[1] ** 2] for r in rows])
    y = np.array([r[3] for r in rows])
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    return coef, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="A,B,C,D")
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--out", default="results/extensions/zakat")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    fn = dict(A=part_a, B=part_b, C=part_c, D=part_d)
    with Pool() as pool:
        for part in a.parts.split(","):
            f = out / f"part_{part}.json"
            if f.exists():
                res[part] = json.loads(f.read_text())
                continue
            print("part", part, flush=True)
            res[part] = fn[part](pool, a.seeds)
            f.write_text(json.dumps(res[part]))
    make_figures(out, res)
    if "A" in res and "B" in res:
        coef, rows = fit_critical(res)
        pred = np.array([coef[0] * r[0] * r[1] * r[2] + coef[1] * r[0] * r[1] ** 2 for r in rows])
        y = np.array([r[3] for r in rows])
        r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        (out / "critical_fit.json").write_text(json.dumps(dict(a=float(coef[0]), b=float(coef[1]), r2=float(r2),
                                                               max_rel_err=float(np.abs(pred / y - 1).max()), n=len(rows))))
        print(f"z_c = {coef[0]:.4f} n_y beta waa + {coef[1]:.4f} n_y beta^2   R2={r2:.4f}")
    print("done")


if __name__ == "__main__":
    main()
