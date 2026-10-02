# Replication of Kato (2022), PLOS ONE 17(9): e0275113

Code: `islamic_econophysics/kato2022.py` (engine), `experiments/replicate_kato2022.py`
(Figs 3-7), `tests/test_kato2022.py`. Outputs: `results/kato2022/` (figures,
`results.json`, `comparison.md`). Run with `python experiments/replicate_kato2022.py`
(~25 s, 4 cores) and `pytest tests`.

The paper releases no code or tables, so "reproduces" is judged against values taken from the
printed figures (`islamic_econophysics/kato2022_reference.py`; extracted from the PDF's
305 dpi embedded images by `experiments/digitize_kato_figures.py`, see "What the PDF does and
does not contain"). Kato shows **one run per point**; this repo averages 20 seeds, so every
Gini comparison is made two ways (`results/kato2022/comparison.md`): a fixed +/-0.05 band, and
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
   `--insolvency skip` (results in `results/kato2022_skip/`).
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

## Bankruptcy options (beyond Kato; L model, `insolvency="bankrupt"`)

At a default the borrower's assets are `X = lam m_i + (1-lam)(1+d)(m_i+m_j)` and the
lender's claim is `C = (1-lam)(1+rho) m_j > X`. Plain liquidation (the default) gives
the lender all of `X`, which reproduces Eq 2's total exactly. Eleven optional knobs, all off
by default (`Params`, off = previous results to 4e-14):

| option | mechanism | wealth |
|---|---|---|
| `recovery=r` | lender settles for `min(r*C, X)`; borrower keeps the rest (haircut / debt relief, cf. Q 2:280) | conserved |
| `bankruptcy_cost=c` | fraction `c` of `X` destroyed before settlement | destroyed |
| `exclusion_period=E` | defaulter cannot borrow for `E` steps (can still lend) | conserved |
| `lender_exclusion_period=E` | defaulter also cannot *lend* for `E` steps | conserved |
| `garnishment_rate=g` | **debt ledger**: lender seizes the borrower's assets as before, and the unpaid part `U` becomes a debt of the borrower to that lender. Whenever the debtor's wealth rises (profit as borrower, interest as lender, transfers received) a fraction `g` of the *gain* goes to the creditor until `U` is repaid. One open debt per borrower; a later lender of an already-indebted borrower just absorbs its loss. Holdings stay >= 0 and the Gini is on holdings, as in the paper. `g=1` is the nearest analogue of literal Eq 2, where `U` stays as negative wealth. Cannot be combined with `recovery`/`bankruptcy_cost`. | conserved (transfers only) |

Ledger extensions (all need `garnishment_rate > 0`, all off by default; each version reproduces
the previous one to roundoff, ~1e-14):

| option | mechanism |
|---|---|
| `multi_debt=True` | debts are an N x N matrix: a borrower can owe several lenders at once, every default adds to what is owed to that lender, and a repayment is split pro rata over the creditors. (Single mode: a borrower with an open debt cannot take a second; the later lender just absorbs the loss.) |
| `seize_receivables=True` | the defaulter's estate includes the claims others owe it. Its lender takes those claims over pro rata up to the unpaid amount `U`; if the lender itself owes the defaulter, that part is set off. This is how a creditor becomes exposed to its debtor's other debtors. (Lenders in this model cannot go negative on their own, so this is the available form of "creditors who go bust in turn".) Needs `multi_debt`. |
| `debt_term=T` | a debt is forgiven (the lender writes it off) `T` steps after the most recent default on that lender, enforced whenever a debt is touched and at every sample. Clocks are per (debtor, lender) pair, so a repeat default restarts the clock. 0 = never. |
| `repayment_order` | how a debtor's repayment is split over its creditors: `pro_rata` (pari passu, default), `fifo` (oldest default first, i.e. seniority by date) or `lifo` (newest first). Needs `multi_debt`. |
| `debt_interest=r` | open debts grow by `(1+r)` per step, compounded lazily when touched and at every sample. Agents are drawn about every N/2 = 500 steps, so `r=1e-4` is roughly 5% per loan event. Without `debt_term` the debts of near-zero agents compound without bound. |
| `receipts_are_income=True` | **funding side**: a creditor that is itself indebted passes a share `g` of each repayment it receives on to its own creditors, cascading down chains of debt (at most 32 hops per income event; the remainder is counted in `cascade_truncated`). |

Diagnostics on `Result`: `recovery_rate` (lender receipts at default plus later repayments /
claims; 0.93 under plain liquidation), `debt_outstanding`, `debt_repaid`, `debt_forgiven`,
`claims_seized`, `neg_net_worth_share` (agents whose holdings + receivables - debts < 0),
`idle_steps` (steps where bans left no eligible pair; resampling is bounded so a ban cannot
hang a run), and two default counts. **`defaults` counts every step on which Eq 2 would
leave the borrower negative, including trivial ones between two agents holding almost nothing;
`material_defaults` counts only those with an unpaid part above 1% of initial mean wealth.**
Use the latter (or `recovery_rate`) for conclusions; the raw count is dominated by near-zero
agents (e.g. 2.4e5 raw vs 1.1e5 material under plain liquidation).

`experiments/bankruptcy_sweep.py` scores each setting against Kato's L curves
(`results/bankruptcy_sweep/sweep.md`, 10 seeds). Gini at t=1e6, rho=0.05 (Kato: 0.98 alone,
0.84 with T at tp=1e5):

| setting | L alone | L-T, tp=1e5 | material defaults | note |
|---|---|---|---|---|
| plain liquidation | 0.94 | 0.92 | 1.1e5 | |
| recovery 0.95 / 0.90 / 0.75 | 0.87 / 0.74 / 0.55 | 0.86 / 0.74 / 0.55 | (raw 2.1e5 / 1.1e5 / 4.0e4) | |
| recovery 0 (lender gets nothing) | 0.49 | 0.49 | 2.6e4 | |
| cost 0.1 / 0.5 | 0.91 / 0.91 | 0.88 / 0.91 | 1.2e4 / 2.8e3 | total wealth -> ~0 |
| exclusion 1e3 / 1e4 | 0.95 / 0.98 | 0.93 / 0.97 | 9.6e4 / 4.1e4 | |
| lender exclusion 1e3 / 1e4 | 0.94 / 0.94 | 0.92 / 0.94 | 1.5e5 / 4.9e5 | |
| ledger (single) g = 0.1 / 0.5 / 1 | 0.94 / 0.96 / 0.97 | 0.93 / 0.95 / 0.96 | 1.1e5 / 9.4e4 / 9.0e4 | 79-83% of agents net-negative |
| ledger multi g=1 | 0.98 | 0.97 | 6.4e4 | debt 29 per agent; 90% net-negative |
| ledger multi+seize g=1 | 0.98 | 0.96 | 7.3e4 | |
| ledger multi, term 1e4 / 1e3 | 0.95 / 0.94 | 0.93 / 0.92 | 1.0e5 / 1.1e5 | debt 0.28 / 0.03 per agent |
| ledger multi, fifo / lifo | 0.97 / 0.97 | 0.96 / 0.95 | 7.7e4 / 6.6e4 | pro rata: 0.98 / 0.97 / 6.4e4 |
| ledger multi, interest 1e-5 (no term) | 0.97 | 0.96 | 7.1e4 | open debt 4.1e4 per agent |
| ledger multi, term 1e4 + interest 1e-4 | 0.95 | 0.93 | 1.0e5 | debt 0.50 per agent; 6.6e4 forgiven |
| ledger multi + receipts (funding side) | 0.98 | 0.97 | 6.1e4 | 307 passed on of 4037 repaid |
| multi + seize + receipts + term 1e4 | 0.95 | 0.93 | 1.0e5 | no chains left to cascade |

What this shows:

- **Partial recovery is a strong equaliser, on its own.** A 10% haircut cuts the exchange-only
  L Gini from 0.94 to 0.74, a 25% haircut to 0.55, with no redistribution at all. The
  defaulter keeps a stake instead of being reset to 0, which is the likely reason (I did
  not isolate it), and the lender recovers less of what it is owed. (Within this model; there
  is no credit-supply side, so nothing here says whether lenders would lend at all under a
  haircut regime.)
- **It does not rescue the Kato fit.** L-T tracks L alone for every r, so no r gives both
  Kato's 0.98 (alone) and 0.84 (with T). The mean deviation over the 45 digitised L points
  is 0.056 for plain liquidation and 0.20-0.26 for any r <= 0.75, i.e. worse.
  Cost 0.1-0.25 (0.048) is within noise of baseline.
- **Bankruptcy cost drives total wealth to ~0** by t=1e6 (there is no production in
  these models, only mean-zero shocks, so deadweight losses accumulate). Gini is scale-free so
  the shares are still meaningful, but read the level of wealth as meaningless there.
- **Exclusion** lowers defaults and pushes the late-time Gini up towards Kato's
  rho=0 curve (0.55 -> 0.69 at E=1e4; L alone 0.98 matches) but overshoots L-T (0.97 vs 0.84).
- **Lender exclusion barely moves inequality** (L alone 0.94 whatever E) but multiplies
  material defaults 4.5x at E=1e4 (1.1e5 -> 4.9e5). My earlier guess was that removing
  defaulters from the lender pool leaves richer lenders whose claims are large relative to
  their borrowers; that is a hypothesis, not tested.
- **The single-debt ledger recovers little because defaulters have no income.** A defaulter is
  reset to ~0 and mostly meets rich lenders, so it rarely gains; open debt sits at about one
  unit per agent, comparable to all wealth in the system (0.85 per agent). 40% (g=0.1) to
  66% (g=1) of the debt created had been repaid by t=1e6. With transfers, a debtor's
  received transfers are partly passed on to creditors, so the ledger works against T.
  The fall in raw defaults at g=1 (2.4e5 -> 1.4e5) is mostly trivial defaults disappearing;
  material defaults fall 1.1e5 -> 0.9e5.
- **Multi-creditor debt overhang.** With every default recorded, open debt grows to ~29 per
  agent, about 35x the wealth in the system, and ~90% of agents have negative net worth,
  while holdings Gini reaches **0.98 for L alone (Kato: 0.98) and 0.75 for L at rho=0
  (was 0.55; Kato 0.95)**. Of all the options this is the one that moves the late-time L
  curves closest to Kato's. L-T becomes 0.97 (Kato 0.84), so the overall mean deviation
  is still 0.077 (baseline 0.056). Caveat: the holdings Gini ignores what the poor owe, so
  on net worth the picture is far worse than 0.98.
- **Seizing receivables has almost no effect on the Gini** (0.98; 2.3e3 claims seized out of
  ~3e4 debt created), because the defaulter's debtors are mostly other near-zero agents.
  It cuts raw defaults (2.9e5 -> 1.1e5), which I read as a change in how many agents sit at
  exactly zero rather than anything economic; material defaults move 6.4e4 -> 7.3e4.
- **Debt expiry works as a bankruptcy discharge.** With term 1e4 open debt falls from 29 to
  0.28 per agent, net-negative agents from 90% to 66% (19% at term 1e3) and the Gini back
  to 0.95 / 0.94. About 95% of all debt ever created is forgiven rather than repaid
  (term 1e4: 2.9e4 forgiven vs 1.4e3 repaid), so in this model "repayment from future
  income" is mostly a fiction and the discharge rule does the work.
- **Creditor priority makes no visible difference to inequality** (fifo 0.97, lifo 0.97,
  pro rata 0.98 for L alone) and only shifts how many defaults are material (7.7e4, 6.6e4,
  6.4e4). Repayments are small relative to debts, so who is paid first rarely matters.
- **Interest on debts only compounds the overhang.** At 1e-5 per step (no term) open debt
  reaches 4.1e4 per agent, about 47 000 times the wealth in the system, with the holdings Gini
  unchanged at 0.97; with a 1e4-step term the same idea is bounded (0.5 per agent) and 6.6e4 of
  debt is forgiven, versus 2.9e4 without interest. Interest on the debts of people with no
  income is a bookkeeping fiction in this model: it raises claims, never recoveries.
- **The funding side (receipts passed down the chain) barely fires.** Only 307 of 4037 units
  repaid (7.6%) were passed on, because few creditors are themselves in debt to anyone who is
  still owed; once debts expire (term 1e4) it is 0. Gini is unchanged (0.98). (A bug found
  while building this: float roundoff in cascades left holdings at about -1e-15, which a
  negative *lender* compounds at (1-lam)*rho per loan, reaching -700 by t=3e5. Repayments are
  now capped at current holdings; the old configurations are unchanged to 1e-14 relative.)
- **Net:** none of the options improves the overall fit to Kato (best mean deviation:
  bankruptcy cost 0.1-0.25, 0.048 vs 0.056, within noise; term 1e3 ties baseline at
  0.057). They are policy levers for the Islamic-finance extensions, not repairs to the
  replication.

Still not built: a credit-supply side (lenders deciding whether to lend given their own
balance sheet, which is what a haircut or discharge regime would really test), secured or
subordinated debt classes other than by date, and debtors' income from anything other than
the model's own exchanges and transfers.

## Corrections to the literature dossier (`overview.md`, which is a PDF despite the extension)

- It says the L and J equations "appear as images in the HTML version"; I checked the PDF
  and the forms above are Kato's. The dossier's reconstruction is correct.
- "Q transfers to the poorest 1/4 to 1/6 of agents" is loose: it is agents holding at
  most `m_max/q`, which is a wealth threshold, not a head-count.
- Paper typos worth knowing: Fig 4's inset says lam=0.5 (the caption and the curves say
  0.25); Fig 3's caption says t up to 1e5 but shows 1e6; the text says Q at 1/1 equals
  "R, L and T" where "R, L and J" is meant. In Fig 4 the legend leader lines for
  "L: rho=0" and "J: dw=0.2" cross (checked at 300 dpi).

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
