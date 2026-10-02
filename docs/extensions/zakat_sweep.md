# Zakat-rate sweep: yard sale + wealth-attained advantage

> **EXTENSION: not part of the Kato (2022) replication.** The paper has no such model. The
> replication status is in `docs/replication/kato2022_replication.md` and `README.md`.

Code: `islamic_econophysics/extensions/yardsale.py` (engine), `experiments/extensions/zakat_sweep.py` (sweep),
`tests/extensions/test_yardsale.py`. Results: `results/extensions/zakat/` (`part_A..D.json`, `critical_fit.json`,
`dose_response.png`, `critical_scaling.png`). Run: `python experiments/extensions/zakat_sweep.py`
(~35 min on 4 cores; parts are cached).

**This is not a replication of Kato (2022), and Kato's paper has no zakat model to replicate.**
He models only a periodic transfer T and a quantile transfer Q. This sweep uses Boghosian's
yard-sale model and the zakat rule from the dossier (2.5% above a nisab per lunar year), so it
cannot be checked against his figures. The L-model uncertainties do not affect it. What *would*
need the author is how his T and Q (and his time step) correspond to zakat; see
`docs/replication/email_to_kato.md`, questions 8-10.

## Model

Each step a random pair trades; the stake is `beta * min(w_i, w_j)` (the poorer agent's wealth
is at risk) and agent i wins with probability `1/2 (1 + waa (w_i - w_j)/(w_i + w_j))`
(waa = 0 is the fair coin; waa > 0 is the wealth-attained advantage). Every lunar year
(`trades_per_agent_year` trades per agent) each agent holding at least `nisab` x mean wealth
pays `zakat_rate` x its wealth, and the pool goes to all agents equally, or only to the agents
below the nisab ("poor"). Everyone starts with wealth 1; total wealth is conserved.
Defaults: N = 2000, beta = 0.1, 20 trades per agent-year, 6000 years, 6 seeds, steady state =
last quarter of the run.

The mapping of `waa`, `zakat_rate` to the continuous-time (zeta, chi) of Boghosian et al.
is **not assumed**; the critical rate is measured.

## What "critical" means here

Boghosian's theory has a sharp transition at chi = zeta in the continuum limit. In the agent
model with finite N the response to the levy is smooth, and the share of the single richest
agent is not a fixed fraction as N grows (it shrinks like a power of N), so "an oligarch
holds half the wealth" has no N-independent threshold. What is N-independent is the **Pareto
exponent alpha of the stationary wealth tail** (Hill estimator on the top 5% of agents): at
alpha <= 1 the mean is carried by the extreme tail, so a few agents hold a finite share
however large N is. I define the critical zakat rate `z_c` as where alpha crosses 1
(log-interpolated). This is my operational definition, not Boghosian's transcritical
bifurcation, though it behaves like one (below).

## Results

**Dose-response** (waa x zakat rate, N=2000, 20 trades per agent-year, beta=0.1;
`dose_response.png`). Gini falls smoothly with the levy from ~0.95-0.99 at 0.05% to ~0.35-0.45
at 10%. At 2.5% per lunar year:

| waa | Gini | top 1% share | top 10% share | bottom 50% share | tail alpha | z_c |
|---|---|---|---|---|---|---|
| 0 | 0.57 | 0.07 | 0.42 | 0.123 | 3.24 | 0.14% |
| 0.02 | 0.59 | 0.08 | 0.44 | 0.113 | 3.06 | 0.21% |
| 0.05 | 0.62 | 0.09 | 0.47 | 0.099 | 2.88 | 0.28% |
| 0.1 | 0.67 | 0.11 | 0.53 | 0.082 | 2.50 | 0.43% |
| 0.2 | 0.73 | 0.14 | 0.63 | 0.060 | 2.03 | 0.72% |

So with these trade parameters 2.5% is 3.5x (waa=0.2) to 17x (waa=0) above the critical
rate: no oligarchy, but still a Gini of 0.6-0.7 (the levy does not make the economy equal).

**Scaling** (`critical_scaling.png`). Over 18 (trades per agent-year n_y, stake beta, waa)
points,

    z_c ~= n_y * beta * (0.0129 * waa + 0.0085 * beta)       (R^2 = 0.999, worst point 18% off)

- **Linear in trade frequency** (z_c = 0.0019, 0.0065, 0.027 at n_y = 5, 20, 80 for waa=0.2):
  a levy per *year* is only meaningful relative to how often agents trade per year.
- The waa term is the advantage's drift per year (n_y * beta * waa); the beta^2 term is the fair
  game's own diffusion, which vanishes in the small-stake limit, as a continuous theory
  requires (it is why z_c(waa=0) is not 0 here). The ratio of the two, ~0.66 beta/waa, says the
  advantage dominates once waa exceeds ~beta.
- **N-independent**: z_c(waa=0.1) = 0.414%, 0.425%, 0.426% at N = 500, 2000, 8000.

**When does 2.5% suffice?** From the fit, the largest trade frequency (trades per agent per
lunar year) for which 2.5% still lies above z_c:

| beta \ waa | 0 | 0.05 | 0.1 | 0.2 | 0.5 | 1 |
|---|---|---|---|---|---|---|
| 0.05 | 1174 | 466 | 291 | 166 | 73 | 37 |
| 0.1 | 294 | 167 | 117 | 73 | 34 | 18 |
| 0.2 | 73 | 53 | 42 | 29 | 15 | 9 |

(Direct check: at beta=0.1, waa=0.2 the measured z_c is 2.7% for 80 trades per agent-year, so
2.5% is *below* critical there.) The answer is therefore conditional on three quantities
that are not observable from this model: trade frequency, stake size and the strength of
the wealth advantage. I do not claim a verdict for any real economy.

**Nisab and recipients** (`part_C.json`, N=2000, 2.5%). The nisab (0-1 x mean wealth) does
not change z_c (0.41-0.45% for every setting) or the tail exponent: it only exempts agents
who were not paying much anyway. Paying the pool only to agents below the nisab (instead of
equally to all) leaves z_c unchanged but moves wealth to the bottom: at waa=0.05 the bottom
half's share rises from 9.9% to 15.6-16.4% (nisab 0.1-0.5 of mean) and the Gini drops from
0.62 to 0.53-0.54; at waa=0.2 from 6.0% to 9.0-10.9% and Gini 0.73 to 0.66-0.68. In this model
who *receives* decides the poor's share; who *pays* (the tail) decides whether oligarchy
forms.

## Caveats

- Toy units. `waa`, `beta` and trades per year are not calibrated to any data; the
  structure (z_c linear in trade frequency and advantage) is the finding, not the numbers.
- z_c is defined by my tail-exponent criterion (Hill, top 5%, 6 seeds, grid interpolation);
  a different tail fraction or horizon shifts the number somewhat. Relaxation takes ~1/z
  years, so rates below ~0.1% are not fully equilibrated at 6000 years.
- Zakat is levied once a year on the whole wealth above the nisab (the majority juristic
  reading); it is not modelled by asset class (stock vs flow, livestock, crops), nor with
  the eight categories of recipients, nor with the interaction of exchange rules (riba
  prohibition, profit-loss sharing) that the L and J models treat.
- The result is about *which levy rate stops tail condensation*, not about welfare.
