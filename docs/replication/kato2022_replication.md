# Replication of Kato (2022), PLOS ONE 17(9): e0275113

Code: `islamic_econophysics/kato2022.py` (engine), `experiments/replication/replicate_kato2022.py`
(Figs 3-7), `tests/replication/test_kato2022.py`. Outputs: `results/replication/kato2022/` (figures,
`results.json`, `comparison.md`). Run with `python experiments/replication/replicate_kato2022.py`
(~25 s, 4 cores) and `pytest tests`.

The paper releases no code or tables, so "reproduces" is judged against values taken from the
printed figures (`islamic_econophysics/kato2022_reference.py`; extracted from the PDF's
305 dpi embedded images by `experiments/replication/digitize_kato_figures.py`, see "What the PDF does and
does not contain"). Kato shows **one run per point**; this repo averages 20 seeds, so every
Gini comparison is made two ways (`results/replication/kato2022/comparison.md`): a fixed +/-0.05 band, and
"within 0.01 + 2.5 sd of my per-run scatter", which is the right yardstick for a single
realisation.

## Equations as implemented (checked against the PDF, not the HTML)

| | |
|---|---|
| R, Eq 1 | `m_i' = lam m_i + eps(1-lam)(m_i+m_j)`, `m_j' = lam m_j + (1-eps)(1-lam)(m_i+m_j)` |
| L, Eq 2 | `m_i' = lam m_i + (1-lam)(m_i - rho m_j + d(m_i+m_j))`, `m_j' = lam m_j + (1-lam)(1+rho) m_j` |
| J, Eq 3 | `m_i' = lam m_i + (1-lam)(1+d) m_i`, same for j |
| T, Eq 4 | `m_i' = (1-xi) m_i + xi * mean_{j!=i} m_j`, all agents at once every `tp` steps |
| Q, Eq 5-6 | `thr = m_max/q`; above thr: `(1-xi) m_i`; at or below: `m_i + xi*S_q/N_q` |
| Gini, Eq 7 | sorted formula |

In L the interest `(1-lam) rho m_j` is charged on the **lender's** wealth, as the text
says; that is what makes a rich lender ruin a poor borrower (see below).

## What the paper leaves open, and what I chose

1. **L model, insolvent borrower (the big one).** Eq 2 literally lets the borrower's
   wealth go negative and, with rho > 0, it runs away (min wealth about -80 at t=1e5,
   -1e18 at 1e6; `insolvency="none"`). Kato says nothing. Rules tried (`Params.insolvency`):

   | rule | what happens | total wealth (rho=.05, t=1e5 / 1e6) |
   |---|---|---|
   | `bankrupt` (**default**) | borrower defaults: reset to 0, **lender absorbs the shortfall** | 0.99 / 0.95 (conserved) |
   | `skip` | trade rejected if borrower would go negative | x1.9 / x1.4e5 (inflates) |
   | `floor`, `cap` | borrower clamped at 0 / pays what it can | inflate |

   `bankrupt` is the explicit-default model: limited liability for the borrower, loss on
   the lender (who lent `(1-lam) m_j`), and the total wealth then evolves exactly as in
   Eq 2, so normalising total wealth is not needed (it would not matter anyway: every model
   here is homogeneous of degree 1 in wealth and the Gini is scale-free). `Result.defaults`
   counts defaults (about 1.6e4 by t=1e5 and 2.4e5 by 1e6 at rho=0.05).
   Fit: `bankrupt` is best for L alone (0.83/0.92/0.94 at 1e5/3e5/1e6 vs Kato's
   0.78/0.92/0.98) but too unequal in L-T at small transfer rates (Fig 6, tp=1e4,
   xi=0.125: 0.81 vs 0.64). `skip` is the reverse, and is selectable with
   `--insolvency skip` (results in `results/replication/kato2022_skip/`).
2. **J model: one delta per pair or one per agent?** Eq 3 uses one symbol, so a shared
   draw is the default (`shared_delta=True`). Gini cannot tell them apart.
3. **Evaluation time.** Fig 3's caption says t up to 1e5 but its panels, Fig 4 and Fig 5
   run to 1e6; Figs 6-7 give no t. 1e6 reproduces them (Fig 6 at xi=0 equals Fig 4 at 1e6).
4. **When Gini is read relative to a transfer.** Fig 5's tp=1e5 curve equals Fig 4's
   exchange-only curve at t=1e5, so g is read *just before* the transfer fires: the peak
   of a saw-tooth. Right after a transfer g is much lower.
5. **"Quantile" in Q is a fraction of the maximum wealth** (threshold `m_max/q`), not a
   population percentile. q=1 means nobody is above the threshold, i.e. no
   redistribution. The x-axis of Fig 7 is q = 1..10 (labelled "1/q").

## Results

With plain liquidation (default): 17 of 25 compared series are within +/-0.05, and 11 of the
18 Gini series are consistent with single-run noise (`max z` <= ~2.5). Everything that
matches is J, R or a redistribution of J/R, plus L at t <= 1e5. What is *not* consistent with
noise is the L model (z = 12-57) and one point of one J-Q curve (z = 3.0).

- **R:** Gini 0.37-0.38. This is not just "about 0.4": for lam=0.25 the stationary law
  is Gamma with shape 1+3lam/(1-lam)=2, whose Gini is exactly 3/8 = 0.375. (That is also
  the Dragulescu-Yakovenko two-earner value, a coincidence of shape-2 gamma = sum of two
  exponentials; it says nothing about wealth.)
- **J:** both widths match Fig 4 to within 0.01 at every t. Gini rises faster for
  dw=0.2 (0.61 vs 0.33 at 1e5). J is a multiplicative random walk with no transfer
  between partners, so this is lognormal spreading, `sigma^2 ~ 2t(1-lam)^2 dw^2 / (3N)`.
- **J-T, R-T:** Figs 5-6 match to <= 0.04 (e.g. J-T dw=0.1: 0.12 at tp=1e4, 0.38 at 1e5;
  dw=0.2: 0.25 and 0.68 vs Kato 0.66). The fall of g with xi and its saturation near
  0.5 are reproduced. The tp dependence is a step of ~0.25-0.5 per decade, so
  "approximately logarithmic" is a fair description of two points but not tested.
- **L-T:** the tp=1e4 plateau matches (0.32 vs 0.31). At tp=1e5 `bankrupt` gives 0.91 at
  t=1e6 vs Kato's 0.84 (`skip`: 0.86), and the xi-dependence is too flat (Fig 6/7: dev up to
  0.17 at xi=0.125, 0.08 in Fig 7).
- **Q:** J-Q matches Fig 7 within single-run noise: every point of both J curves is within
  1.5 sd of my ensemble except one, J dw=0.1 at q=10 (Kato 0.471, mine 0.40 +/- 0.02, z=3.0).
  An earlier version of these notes called the U-shapes at large q a mismatch; with the
  digitised values and per-run scatter they are not, apart from that one point. The
  fixed +/-0.05 band would flag J dw=0.2 (0.056) although its max z is only 1.5.

### Not reproduced

- **L at rho=0, t >= 3e5** (Fig 4). Kato: 0.75 at 3e5, 0.95 at 1e6. Mine (any rule): 0.55.
  Without interest the only force on a poor borrower is the zero-mean shock
  delta*(m_i+m_j), which is dominated by the rich lender's wealth, so a reflecting floor at 0
  (default and lender write-off) caps the Gini. Kato's values need either negative wealth or
  something he did not describe. Up to t=1e5 it matches (0.27, 0.45 vs 0.49).
- **L-T / L-Q shapes** under `bankrupt` (above), and the L Gini at xi=0 / q=1 reads 0.94
  rather than 0.98.
- **Fig 3's L panel, t=1e6** (a quantitative target I had not used): Kato has ~885 of 1000
  agents in the first 0.1-wide wealth bin and almost none in the next; J matches (dw=0.1: 408
  vs ~390; dw=0.2: 906 vs ~910) but plain-liquidation L has 741 (and a tail of 48, 29, 20 in
  the next bins). Raw wealth is plotted, so this also pins the level of total wealth:
  Kato's t=1e3/1e4 panels centre on m=1, i.e. the L model's wealth does not blow up or
  collapse. Options that land nearer: single-debt ledger g=1 (866), exclusion 1e4 (912),
  multi-debt ledger g=1 (936), multi g=0.5 (827). Bankruptcy cost drives all 1000 below 0.1
  (wealth -> 0). None of them also gets L-T right (Fig 5: 0.84).

## Beyond the paper: bankruptcy options and debt ledger (moved)

The L model's missing insolvency rule led me to build a set of bankruptcy options (partial
recovery, bankruptcy cost, borrowing and lending bans, a debt ledger with priority, interest, expiry,
receivable seizure and a funding cascade). **None of that is in Kato (2022)** and none of it is used
by the replication (the replication script asserts that every such option is off). It moved to
`docs/extensions/l_bankruptcy.md`; the results are in `results/extensions/l_bankruptcy/`.

## Corrections to the literature dossier (`overview.md`, which is a PDF despite the extension)

- It says the L and J equations "appear as images in the HTML version"; I checked the PDF
  and the forms above are Kato's. The dossier's reconstruction is correct.
- "Q transfers to the poorest 1/4 to 1/6 of agents" is loose: it is agents holding at
  most `m_max/q`, which is a wealth threshold, not a head-count.
- Paper typos worth knowing: Fig 4's inset says lam=0.5 (the caption and the curves say
  0.25); Fig 3's caption says t up to 1e5 but shows 1e6; the text says Q at 1/1 equals
  "R, L and T" where "R, L and J" is meant. In Fig 4 the legend leader lines for
  "L: rho=0" and "J: dw=0.2" cross (checked at 300 dpi).

## The search for Kato's missing L-model rule

Everything except the L model reproduces within single-run noise. For L I tried, in order,
to find a natural rule that fits *all* the digitised L targets (Fig 3 first bin, Fig 4 rho=.05
and rho=0, Fig 5-7 L-T and L-Q). Summary of what failed and what it taught:

| hypothesis | outcome |
|---|---|
| Eq 2 literal, Gini on raw balances | Gini > 1 or negative (wealth ~ -1e18 by t=1e6); not what is plotted |
| Eq 2 literal, Gini on clipped / abs / positives-only / shifted / squared balances | 0.50-0.76 at t>=1e5 vs 0.98; ~500 agents in the first bin vs 885 |
| interest on the borrower's own wealth (what I first thought Fig 1C showed) | **a misreading**: at 305 dpi Fig 1C shows -rho(1-lam) m_j, same as Eq 2 |
| borrower = richer / poorer agent instead of random | Gini 0.07 (richer) or far too fast (poorer) |
| floor, skip, cap, bankrupt, reflect (abs), redraw delta, swap roles | all match through t=1e5; none reaches 0.98 / 0.95 at t=1e6 (0.80-0.96 / 0.50-0.58) |
| shock proportional to the borrower's *own* wealth (loss base 2 m_i) | fits rho=0 well (0.455, 0.70, 0.935 vs 0.493, 0.756, 0.951) but with interest runs too fast (0.90 at 1e5 vs 0.78) or too slow, and L-T lands at 0.96 |
| mixtures of the two interest bases x two loss bases (25 combinations) | best worst-case deviation 0.116 |

What it taught:

- **Eq 2 is right** at least up to the point where balances would turn negative: the literal
  rho=0 run tracks Kato's curve to t=1e5 (0.27, 0.51 vs 0.276, 0.493).
- **Kato's late rho=0 growth looks multiplicative, not additive.** A log-normal spread with
  sigma^2 ~ 7.5 gives Gini 0.947 at t=1e6 (Kato 0.951). In Eq 2 a poor borrower facing a rich
  lender takes a shock proportional to the *lender's* wealth, which any non-negativity rule
  turns into mean reversion at zero, capping the Gini near 0.55. A shock proportional to the
  borrower's own wealth reproduces the rho=0 curve but not the rest.
- **The best partial fits** (`skip`, `redraw`, `swap`, all selectable with `--insolvency`)
  pin L-T at tp=1e5 to within ~0.02 (swap 0.77/0.82/0.82 vs 0.772/0.838/0.842) but not the
  other L figures. Excluding the no-transfer endpoints, no rule gets Fig 6/7 L within noise
  (max z >= 3.8; `skip` is nearest on Fig 7, 0.022).
- **No rule found fits both** L alone at late times (0.98) **and** L with transfers (0.84)
  **and** Fig 3's ~885 near-zero agents. I have not found, and cannot derive from the paper,
  the rule Kato's code used. A fit with enough free knobs would match, but that is curve
  fitting, not replication, so I have not done it.

Questions whose answers would settle it (for the author, kato.takeshi.3u@kyoto-u.ac.jp, or the
PLOS peer-review file): (1) what happens in the code when Eq 2 would make the borrower's wealth
negative; (2) is the Gini computed on raw balances; (3) is total wealth rescaled during runs;
(4) are Figs 3-7 single runs or averages, and with what seeds; (5) is the J profit/loss
draw shared by the pair.

## What the PDF does and does not contain

Asked whether anything in the PDF is out of my sight and blocking a perfect replication: no.

- **Text and equations.** All of it read; Eqs 1-4 and Fig 1C/D checked visually on the rendered
  pages (pdftotext drops minus signs, so the text layer alone is not enough). The upload is
  byte-identical to `journal.pone.0275113.pdf` in the repo and to arXiv 2206.05443v2.
- **Figures.** They are raster images (305 dpi) embedded in the PDF. I first read values off a
  110 dpi render by eye; I have since extracted the images and digitised every legible
  marker programmatically (y axis calibrated on the tick labels, ~0.002 Gini/pixel). The
  by-eye readings were within 0.010 of the digitised ones (mean 0.004), so reading was not the
  limit. Markers hidden behind other markers (e.g. L at rho=0.05 for t<=1e4 and at 1e6, L-Q at
  q=1) and hollow markers (the tp=1e4 series) cannot be recovered exactly and stay by eye.
- **Run-to-run scatter is visible in the paper itself:** in Fig 6 the xi=0 point of J dw=0.1
  is ~0.83 in the tp=1e5 series and ~0.80 in the tp=1e4 series although the configurations
  are identical (no transfer happens at xi=0). That ~0.03 is Kato's single-run noise; no
  replication can pin his curves tighter than that.
- **arXiv v1 (June 2022, pre-review, 35 pp)** read: same equations, same parameters, nothing
  on how negative wealth, seeds or averaging were handled. (v1's Fig 4 caption says lam=0.25,
  the published figure's inset says 0.5; the curves match 0.25.)
- **Not in the PDF at all:** code, random seeds, number of runs, the treatment of a borrower
  whose Eq 2 balance goes negative, whether J uses one delta per pair, and when exactly Gini
  is read relative to a transfer. These, not missing pixels, are what stop a perfect
  replication. The PLOS peer-review history and author responses are published with the
  article and might settle some of them; I have not been able to read them yet.

## Suggested next steps

1. Bankruptcy options, the extended ledger, priority, debt interest and the funding cascade are
   all in (see above); a credit-supply side is the open extension.
2. Then move to the dossier's Phase 1b: the zakat-rate sweep with a WAA term.
