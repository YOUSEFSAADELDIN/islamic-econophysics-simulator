# Replication vs. digitised Kato (2022) figures (20 seeds/point)

L-model insolvency rule: `bankrupt`.

Tolerance +/-0.05 on the Gini index (figure reading error + Monte-Carlo noise).

| series | max abs dev | mean abs dev | within tol |
|---|---|---|---|
| Fig4 R | 0.016 | 0.010 | yes |
| Fig4 L rho=0.05 | 0.051 | 0.021 | **NO** |
| Fig4 L rho=0 | 0.403 | 0.097 | **NO** |
| Fig4 J dw=0.2 | 0.009 | 0.004 | yes |
| Fig4 J dw=0.1 | 0.009 | 0.004 | yes |
| Fig5 R tp=1e4 plateau | 0.007 | 0.007 | yes |
| Fig5 L tp=1e4 plateau | 0.012 | 0.012 | yes |
| Fig5 J0.2 tp=1e4 plateau | 0.004 | 0.004 | yes |
| Fig5 J0.1 tp=1e4 plateau | 0.005 | 0.005 | yes |
| Fig5 R tp=100000 | 0.016 | 0.010 | yes |
| Fig5 L tp=100000 | 0.078 | 0.036 | **NO** |
| Fig5 J0.2 tp=100000 | 0.029 | 0.013 | yes |
| Fig5 J0.1 tp=100000 | 0.007 | 0.003 | yes |
| Fig6 L tp=10000 | 0.173 | 0.061 | **NO** |
| Fig6 J0.2 tp=10000 | 0.011 | 0.005 | yes |
| Fig6 J0.1 tp=10000 | 0.013 | 0.006 | yes |
| Fig6 L tp=100000 | 0.077 | 0.057 | **NO** |
| Fig6 J0.2 tp=100000 | 0.038 | 0.017 | yes |
| Fig6 J0.1 tp=100000 | 0.015 | 0.008 | yes |
| Fig7 L | 0.076 | 0.063 | **NO** |
| Fig7 J0.2 | 0.057 | 0.023 | **NO** |
| Fig7 J0.1 | 0.073 | 0.026 | **NO** |
