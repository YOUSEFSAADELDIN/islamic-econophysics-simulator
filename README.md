# Islamic econophysics simulator

Two clearly separated things live in this repo:

1. **`replication/`: a re-implementation of Kato (2022)**, *Islamic and capitalist economies:
   Comparison using econophysics models of wealth exchange and redistribution*, PLOS ONE 17(9):
   e0275113 (`journal.pone.0275113.pdf`). Goal: reproduce the paper's results. **Not finished.**
2. **`extensions/`: everything the paper does not contain.** A yard-sale economy with
   wealth-attained advantage and a zakat levy (incl. the eight recipient categories), and a set of
   L-model bankruptcy / debt-ledger options. None of it is checked against Kato's figures, because
   there is nothing in the paper to check it against.

## Status of the Kato replication

Kato released no code and no tables; targets are his figures, digitised from the 305 dpi images in
the PDF. His points are single runs, so agreement is judged against single-run scatter (see
`results/replication/kato2022/comparison.md`).

| Part of the paper | Status |
|---|---|
| R model (Eq 1), Gini (Eq 7) | reproduced (Gini 0.375, the exact Gamma(2) value) |
| J model (Eq 3), Figs 3-4 | reproduced (within 0.01 of every Fig 4 point; Fig 3 bar heights match) |
| T and Q redistribution (Eqs 4-6) with R and J, Figs 5-7 | reproduced within single-run noise, except one point (J-Q, dw = 0.1, q = 10, 3 sigma) |
| **L model (Eq 2), exchange only** | matches up to t of about 1e5 (e.g. rho = 0: 0.27, 0.51 vs 0.276, 0.493). **Not reproduced at late times** (0.98 / 0.951 at t = 1e6 vs 0.80-0.96 / 0.50-0.58 for the simple rules) |
| **L with redistribution (L-T, L-Q), Figs 5-7; Fig 3 L panel** | **not reproduced** (several rules fit pieces, none fits all) |

So **the replication is not complete**: every failing item involves the L model. The paper does not say
what happens when Eq 2 gives a borrower a negative balance, and under every rule I tried the
late-time behaviour differs from his. The unknown is a detail of his code, not something visible in
the PDF. The questions that would settle it are in `docs/replication/email_to_kato.md`; the search
for the missing rule is in `docs/replication/kato2022_replication.md`.

## Layout

```
islamic_econophysics/
  kato2022.py             replication engine (R, L, J, T, Q, Gini). EXTENSION_PARAMS lists the options that
                          go beyond the paper (off by default; the replication runs with them off)
  kato2022_reference.py   figure values digitised from the PDF
  extensions/yardsale.py  yard-sale + WAA + zakat (+ eight categories)       [extension]
experiments/
  replication/            replicate_kato2022.py (Figs 3-7 vs the paper), digitize_kato_figures.py
  extensions/             bankruptcy_sweep.py, zakat_sweep.py, zakat_categories.py
tests/replication/  tests/extensions/
docs/
  replication/            kato2022_replication.md (method, results, open points), email_to_kato.md (draft)
  extensions/             l_bankruptcy.md, zakat_sweep.md, zakat_categories.md
results/replication/  results/extensions/
```

## Running

```
pip install -r requirements.txt
pytest tests                                               # ~20 s
python experiments/replication/replicate_kato2022.py       # Figs 3-7 vs Kato, ~25 s
python experiments/extensions/zakat_sweep.py               # extension, ~35 min
```

`overview.md` is the literature dossier (it is a PDF despite the extension).
