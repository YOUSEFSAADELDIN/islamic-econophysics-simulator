# Replication of Kato (2022), PLOS ONE 17(9): e0275113

Code: `islamic_econophysics/kato2022.py` (engine), `experiments/replicate_kato2022.py`
(Figs 3-7), `tests/test_kato2022.py`. Outputs: `results/kato2022/` (figures,
`results.json`, `comparison.md`). Run with `python experiments/replicate_kato2022.py`
(~25 s, 4 cores) and `pytest tests`.

The paper releases no code or tables, so "reproduces" means "within +/-0.05 Gini of
values read off the printed figures" (`islamic_econophysics/kato2022_reference.py`).
The tolerance covers eyeballing error plus Monte-Carlo noise; Kato shows one run per
point, this repo averages 20.

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
   wealth go negative and, with rho > 0, it runs away (min wealth about -80 at t=1e5, -1e18 at 1e6;
   `insolvency="none"`). Kato says nothing. I tried three rules; fit to Figs 4-7:
   `floor` (clamp at 0) and `cap` over-shoot at t=1e5 (0.85 vs 0.78); `skip` (reject
   the trade) matches best (0.79 at 1e5, 0.85 at 3e5) and is the default. **Caveat:**
   `skip` keeps gains but rejects ruinous losses, so total wealth drifts up (x1.9 at
   t=1e5, x1.4e5 at 1e6 for rho=0.05). Gini is scale-free so Figs 5-7 are unaffected
   in kind, but it is a modelling artefact, not Kato's mechanism.
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

15 of 22 compared series are within +/-0.05 (`results/kato2022/comparison.md`).

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
- **L-T:** matches Fig 5 for t <= 1e6 (0.86 vs 0.84 at tp=1e5; 0.32 vs 0.31 at 1e4).
- **Q:** same qualitative picture as Fig 7 (J-Q with dw=0.1 reaches the J-T level at
  q=5-7), but two features differ: mine does not turn up at q=10 (0.40 vs Kato's 0.47)
  and J-Q dw=0.2 falls monotonically (0.75 at q=10) where Kato's has a shallow minimum
  at q~6-7. Max deviation 0.06-0.07.

### Not reproduced

- **L alone at late times** (Fig 4, t >= 3e5). Kato: Gini 0.95-0.98 at 1e6 for both rho.
  Mine with `skip`: 0.88 (rho=0.05) and 0.54 (rho=0), the latter because poor borrowers'
  trades are rejected and the system freezes. With `none` the rho=0 run tracks Kato until
  1e5 (0.27, 0.50) and then goes negative. Everything up to t=1e5, which covers the
  interest-vs-risk comparison (L at rho=0 is clearly worse than J at the same width),
  does reproduce. The same gap makes Fig 6/7 L points at xi=0 / q=1 read 0.88, not 0.98.
- Fig 3's L panel at t=1e6 (Kato: ~900 of 1000 agents near 0). Mine are spread over
  five orders of magnitude, because of the wealth inflation in point 1.

## Corrections to the literature dossier (`overview.md`, which is a PDF despite the extension)

- It says the L and J equations "appear as images in the HTML version"; I checked the PDF
  and the forms above are Kato's. The dossier's reconstruction is correct.
- "Q transfers to the poorest 1/4 to 1/6 of agents" is loose: it is agents holding at
  most `m_max/q`, which is a wealth threshold, not a head-count.
- Paper typos worth knowing: Fig 4's inset says lam=0.5 (the caption and the curves say
  0.25); Fig 3's caption says t up to 1e5 but shows 1e6; the text says Q at 1/1 equals
  "R, L and T" where "R, L and J" is meant. In Fig 4 the legend leader lines for
  "L: rho=0" and "J: dw=0.2" cross (checked at 300 dpi).

## Suggested next steps

1. Decide how to treat insolvency in L, since it carries the exchange-only headline
   ("interest-bearing loans -> Gini 1"). Options: keep `skip` and report only t <= 1e5;
   model bankruptcy explicitly (debt written off against the lender); or add wealth
   normalisation as Kato suggests.
2. Then move to the dossier's Phase 1b: the zakat-rate sweep with a WAA term.
