# Zakat pooled over the eight recipient categories (Q 9:60)

Code: `islamic_econophysics/yardsale.py` (`recipients="eight"`, `CATEGORIES`, `PROFILES`),
`experiments/zakat_categories.py`, `tests/test_yardsale.py`. Results: `results/zakat_categories/`
(`tables.md`, `part_A..E.json`, `league_table.png`, `poverty_vs_rate.png`).
Run: `python experiments/zakat_categories.py` (~10 min, 4 cores).

Same independent yard-sale economy as `docs/zakat_sweep.md` (nothing here depends on Kato 2022).
Base case: N=2000, beta=0.1, waa=0.1, 20 trades per agent-year, 6000 lunar years, 6 seeds,
zakat 2.5% on wealth above a nisab of 0.5 x mean, poverty line 0.2 x mean.

## How each category is modelled

The levy is unchanged (2.5% of the wealth of every agent above the nisab). The pool is split by
`category_weights`; eligibility is decided on wealth *after the levy and before any grant*
(a first version let someone lifted out of the poor band by the *fuqara* grant be paid again as
*masakin*; a unit test caught it). Shares that cannot be spent (nobody eligible) are returned to
everyone equally and reported as "unspent".

| # | category | model | status |
|---|---|---|---|
| 1 | *fuqara*, the poor | wealth < poverty line; paid pro rata to the gap to the nisab, never beyond it | direct |
| 2 | *masakin*, the needy | poverty line <= wealth < nisab; same rule | direct |
| 3 | *amilin*, administrators | wages: equal split over 1% of agents drawn at random each year | direct (wage as pure transfer; no real overhead) |
| 4 | *mu'allafat al-qulub* | equal split over 5% of agents drawn at random from those below mean wealth | **interpretive** |
| 5 | *riqab*, freeing captives | agents below the poverty line for >= 5 consecutive years, equal split | **interpretive** (poverty-trap release) |
| 6 | *gharimin*, debtors | agents who lost > 50% of last year's wealth and are now below the nisab; paid pro rata to the loss, never above it | **interpretive**: a yard-sale economy has no loan book; this is distress relief. Ties to Kato's L-model ledger only in spirit |
| 7 | *fi sabilillah* | equal per-capita grant to all (public good); also takes any unspent share | direct |
| 8 | *ibn al-sabil*, wayfarer | 2% of agents drawn at random, any wealth, equal split | **interpretive** |

Allocation profiles (`PROFILES`) are **stylised, not fiqh rulings**: `all_public` (the earlier
`recipients="all"`), `poor_only` (50/50 fuqara/masakin), `eight_equal` (1/8 each), `seven_equal`
(no *mu'allafat*), `needy_first` (35/35 poor and needy, 10 administrators, 10 debtors, 5 riqab, 5
public). How the schools actually weight the categories (the dossier notes the classical text
fixes the eight recipients and leaves proportions to the jurists) should be taken from
al-Mughni and al-Qaradawi's *Fiqh al-Zakah* before any of these labels are used in the app.

## Results (z = 2.5%, nisab 0.5; `results/zakat_categories/tables.md`)

Poverty headcount = share of agents below 0.2 x mean wealth, measured at the end of the trading
year **before zakat arrives** (the worst point of the cycle). Measured just after, a rule that
tops everyone up scores 0 by construction (my first headline did exactly that and was wrong).

| allocation | Gini | poverty before zakat | bottom-50% share | pool reaching poorest 20% |
|---|---|---|---|---|
| no zakat | 1.00 | 99% | 0.0% | - |
| all_public (equal per capita) | 0.66 | 36% | 8.4% | 20% |
| eight_equal / seven_equal | 0.62 | 28% | 11.3% | 36-38% |
| needy_first | 0.59 | 22% | 13.2% | 46% |
| poor_only | 0.57 | 19% | 14.4% | 52% |
| only fuqara | 0.56 | 16% | 14.9% | 100% |
| only riqab | 0.61 | 27% | 11.7% | 93% |
| only gharimin | 0.63 | 32% | 9.9% | 32% (57% unspent, returned) |
| only masakin | 0.64 | 35% | 9.1% | 4% (22% unspent) |
| only amilin | 0.80 | 62% | 0.1% | 20% |
| only ibn al-sabil | 0.75 | 53% | 1.0% | 20% |

- **Targeting is what matters, not the levy.** The same 2.5% pool cuts poverty from 36% (equal
  grants) to 19% (poor and needy only) to 16% (all to the poor). Splitting it equally over the
  eight categories (28%) leaves poverty about 1.5x the poor-only value; only about a third of
  the pool reaches the poorest fifth.
- **Random-recipient categories can leave inequality worse than equal grants.** Paying the whole
  pool to a handful of random agents (*amilin* as wages, or the wayfarer) gives Gini 0.75-0.80
  and poverty 53-62% (equal grants: 0.66 and 36%): a lottery on a skewed distribution creates new winners without lifting
  the bottom. (A realistic *amilin* share is capped, below.)
- **Administrative overhead costs the poor about 3.5-4.5 points of poverty per eighth** diverted
  (amilin weight 0 / 1/8 / 1/4 / 0.4: poverty 18.9 / 22.5 / 27.0 / 33.3%; pool reaching the
  poorest 20%: 52 / 43 / 34 / 28%). This treats wages as transfers to random agents, so it is a
  *lower* bound on the harm of real overhead, which is partly lost rather than re-spent.
- **The categories that need a funded definition of need (fuqara, masakin) do nearly all the
  work.** Single-category runs rank: fuqara (16%) < riqab (27%) < gharimin (32%) < masakin
  (35%, but only because nobody below the poverty line gets paid) < public (36%) <
  muallafat (36%) < wayfarer (53%) < amilin (62%). The interpretive categories 4-8 add little
  or nothing beyond what a good poverty definition already provides.
- **Zakat rate needed for a given poverty target** (log-interpolated from 0.5/1/2.5/5/10%):

| poverty before zakat falls to | all_public | poor_only | needy_first | eight_equal |
|---|---|---|---|---|
| 20% | 4.8% | 2.4% | 2.6% | 3.4% |
| 10% | 7.7% | 3.3% | 3.6% | 5.7% |
| 5% | > 10% | 4.5% | 4.9% | 7.8% |

  So at the classical 2.5% only the targeted profiles get near 20%; a 10% target needs about
  3.3-5.7% in this economy. These numbers depend on the toy parameters (waa, beta, trade
  frequency) as in the earlier doc; the *ranking* of the profiles does not.
- **The critical rate is unaffected by the split** (tail exponent = 1 at 0.43% / 0.44% / 0.45% for
  poor_only / all_public / eight_equal), consistent with the earlier finding: who pays decides
  whether an oligarchy forms; who receives decides how the poor fare.
- **Nisab** (poverty line fixed): 0.3 / 0.5 / 1.0 x mean gives poverty 17.6 / 18.9 / 27.2% for
  poor_only and 26.6 / 27.9 / 32.2% for eight_equal. A high nisab shrinks the payer base (and
  widens the "needy" band that must be filled from a smaller pool).

## Caveats

- Categories 4, 5, 6 and 8 are my modelling choices (marked above); changing them changes
  their individual rows but not the headline that targeting dominates, because they are
  small shares of the equal split.
- Poverty is relative (a fraction of mean wealth); in a yard-sale economy with total wealth fixed
  there is no absolute poverty. Wages to *amilin* are re-spent as transfers, so overhead is
  understated.
- The *gharimin* rule is a distress proxy; to tie it to real debt, run it on the L model's
  ledger (open debts of agents), which needs the L model's insolvency rule, still unresolved.
