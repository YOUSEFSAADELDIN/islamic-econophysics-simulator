# L-model bankruptcy options and debt ledger

> **EXTENSION: not part of the Kato (2022) replication.** The paper has no such model. The
> replication status is in `docs/replication/kato2022_replication.md` and `README.md`.

Built while looking for the rule Kato's L model uses when a borrower's balance would go negative.
All options are off by default (`kato2022.EXTENSION_PARAMS`); the engine is `islamic_econophysics/kato2022.py`, the sweep
`experiments/extensions/bankruptcy_sweep.py`.

## Options

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

`experiments/extensions/bankruptcy_sweep.py` scores each setting against Kato's L curves
(`results/extensions/l_bankruptcy/sweep.md`, 10 seeds). Gini at t=1e6, rho=0.05 (Kato: 0.98 alone,
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

