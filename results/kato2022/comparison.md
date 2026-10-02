# Replication vs. digitised Kato (2022) figures (20 seeds/point)

L-model insolvency rule: `bankrupt`.

Two criteria. **Band**: max |sim - Kato| <= 0.05 Gini. **Single-run noise**: Kato's points are single runs, so each point must lie within 0.01 + 2.5 sd of my ensemble, where sd is the per-run scatter over seeds; `max z` is the largest |sim - Kato| / sd (z > 3 is not noise).

| series | max abs dev | mean abs dev | within band | max z | consistent with single-run noise |
|---|---|---|---|---|---|
| Fig4 R | 0.016 | 0.010 | yes | 2.6 | yes |
| Fig4 L rho=0.05 | 0.050 | 0.020 | yes | 11.8 | **NO** |
| Fig4 L rho=0 | 0.404 | 0.098 | **NO** | 56.6 | **NO** |
| Fig4 J dw=0.2 | 0.012 | 0.005 | yes | 1.5 | yes |
| Fig4 J dw=0.1 | 0.015 | 0.005 | yes | 3.9 | yes |
| Fig3 agents in first bin, J dw=0.1 (agents) | 22.950 | 22.950 | yes | - | - |
| Fig3 agents in first bin, J dw=0.2 (agents) | 0.550 | 0.550 | yes | - | - |
| Fig3 agents in first bin, L rho=0.05 (agents) | 151.350 | 151.350 | **NO** | - | - |
| Fig5 R tp=1e4 plateau | 0.007 | 0.007 | yes | - | - |
| Fig5 L tp=1e4 plateau | 0.012 | 0.012 | yes | - | - |
| Fig5 J0.2 tp=1e4 plateau | 0.004 | 0.004 | yes | - | - |
| Fig5 J0.1 tp=1e4 plateau | 0.005 | 0.005 | yes | - | - |
| Fig5 R tp=100000 | 0.016 | 0.010 | yes | 2.6 | yes |
| Fig5 L tp=100000 | 0.073 | 0.035 | **NO** | 19.2 | **NO** |
| Fig5 J0.2 tp=100000 | 0.031 | 0.013 | yes | 1.9 | yes |
| Fig5 J0.1 tp=100000 | 0.010 | 0.004 | yes | 7.2 | yes |
| Fig6 L tp=10000 | 0.173 | 0.061 | **NO** | 44.3 | **NO** |
| Fig6 J0.2 tp=10000 | 0.011 | 0.005 | yes | 1.3 | yes |
| Fig6 J0.1 tp=10000 | 0.013 | 0.006 | yes | 3.5 | yes |
| Fig6 L tp=100000 | 0.073 | 0.056 | **NO** | 19.6 | **NO** |
| Fig6 J0.2 tp=100000 | 0.036 | 0.017 | yes | 1.0 | yes |
| Fig6 J0.1 tp=100000 | 0.010 | 0.006 | yes | 1.4 | yes |
| Fig7 L | 0.073 | 0.062 | **NO** | 16.4 | **NO** |
| Fig7 J0.2 | 0.056 | 0.023 | **NO** | 1.5 | yes |
| Fig7 J0.1 | 0.074 | 0.023 | **NO** | 3.0 | **NO** |
