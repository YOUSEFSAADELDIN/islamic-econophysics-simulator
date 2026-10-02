# Draft email to Dr. Takeshi Kato (not sent)

**To:** kato.takeshi.3u@kyoto-u.ac.jp (corresponding author, PLOS ONE 17(9): e0275113)
**Subject:** Replicating your Islamic/capitalist econophysics models – a few methodology questions

Dear Dr. Kato,

I am building an open-source educational simulator of wealth exchange and redistribution
under Islamic and capitalist rules, and I started by re-implementing the models in your 2022
PLOS ONE paper from the published equations (no code was released, and the data statement
says everything is in the paper). I would like to credit your work properly and, if you are
willing, to reproduce it exactly. I am writing once with all of my methodology questions so as
to take as little of your time as possible; one-line answers are plenty.

**What reproduces.** With N = 1000, λ = 0.25 and your parameters, the R, J, T and Q models
match your figures to within the scatter of a single run: R at Gini 0.375 (the exact value of
the Gamma(2) stationary law), J at δw = 0.1 and 0.2 within 0.01 at every point of Fig 4,
J-T and R-T (Figs 5-6) and J-Q (Fig 7) within single-run noise (one point, J-Q δw = 0.1 at q = 10, is 3 standard deviations off: 0.47 vs my 0.40), and your Fig 3 histograms for
J (first-bin counts 408 vs about 390, and 906 vs about 910 at t = 10^6). I digitised your
figures from the embedded images to compare at about 0.002 Gini resolution.

**What does not: the L model at late times.** Equation 2 as printed reproduces your curves
up to t of about 10^5 (e.g. ρ = 0: 0.27 and 0.51 at t = 3x10^4 and 10^5 against your 0.276 and
0.493). After that, the borrower's balance in Eq 2 can become negative (a poor borrower facing
a rich lender takes a shock proportional to the lender's wealth, and with ρ > 0 pays interest
on the lender's wealth), and I cannot tell what your code does there. Taking Eq 2 literally
the wealth diverges (about -10^18 by t = 10^6). The simple non-negativity rules I tried (floor at zero, reject the trade, lender absorbs
the shortfall, redraw δ, swap roles, partial recovery) all match up to 10^5 and then stall:
Gini 0.80-0.96 at t = 10^6 for ρ = 0.05 (you: 0.98) and 0.50-0.58 for ρ = 0 (you: 0.951), and
Fig 3B has about 741 of 1000 agents below 0.1 (you: about 885). A debt-ledger variant (the
unpaid balance stays owed and is garnished from the debtor's later gains) reaches 0.98 and
0.75 but then gives 0.97 for L-T at tp = 10^5, where you have 0.84. I could not find a rule
that fits all of Figs 3-7 for L together.

**Questions**

1. *L model, insolvency.* When Eq 2 would give the borrower a negative wealth, what does the
   code do (allow it, floor at zero, reject, something else)? Is the Gini (Eq 7) computed on
   the raw balances?
2. *Rescaling.* Is total wealth rescaled during the L and J runs (you mention it as an option
   in the Results)? If so, how and how often?
3. *J model.* Is the profit/loss rate δ in Eq 3 one draw shared by the pair i, j, or one draw
   per agent?
4. *Runs.* Are the points in Figs 3-7 single runs or averages, and with how many seeds? (In
   Fig 6 the ξ = 0 point of J-T δw = 0.1 is about 0.83 for tp = 10^5 and about 0.80 for
   tp = 10^4 although no transfer happens at ξ = 0, which suggests single runs.)
5. *Timing.* In Figs 5-7, is the Gini read just before the transfer at multiples of tp
   (your tp = 10^5 curve equals the exchange-only curve at t = 10^5)? At what time are Figs 6
   and 7 evaluated (10^6?)?
6. *Q model.* Is the threshold m_MAX/q a fraction of the maximum wealth (not a population
   quantile), and is "1/1 quantile equals the R, L and T models" meant to read "R, L and J"?
7. *Small inconsistencies.* The inset of Fig 4 says λ = 0.5 (your curves and caption fit
   0.25); the caption of Fig 3 says t = 10^3, 10^4, 10^5 but the panels are labelled 10^6.
8. *Zakat.* Your abstract names waqf, sadaqah and zakat as the Islamic redistribution, but the
   models are T (which you equate with inheritance tax and waqf) and Q. Do you regard T/Q as
   representing zakat (2.5% per lunar year of wealth above a nisab, paid to the poor, the
   needy and the other classes of Qur'an 9:60), or was zakat deliberately left outside the
   model? Is there a zakat-specific model, or are there data you would use for it?
9. *Time unit.* What does one time step stand for? You equate t = 10^4-10^5 with a millennium
   ("12 months x 365 days x 1000 years"), which suggests a step is roughly a month or a day.
   With that mapping, which (xi, tp) would you use for an annual 2.5% levy?
10. *Recipients.* Is Q (transfer to agents below m_MAX/q) meant as "the poor" only, and do you
    have a view on how the other recipient classes would enter such a model (for example
    administrators, debtors or public spending)?
11. *Code.* Would you be willing to share the simulation code (Mathematica?) or the random seeds
   for a few runs? Even a single L-model run with the intermediate wealth vector at t = 10^6
   would settle it.

Everything is public (or I can attach it): the equations as implemented, the digitised
figure values, a scored comparison against each figure, and the L-model rule search:
[repo link]. I am glad to share the full write-up (`docs/kato2022_replication.md`).

With thanks for the paper and for any help,

[Name]
[Affiliation / contact]
