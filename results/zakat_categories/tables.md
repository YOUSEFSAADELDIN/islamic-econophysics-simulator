## A. One category at a time, and stylised profiles (z=2.5%, nisab 0.5, N=2000, waa=0.1)

| allocation | Gini | poverty before zakat / after (<0.2 mean) | bottom-50% share | top-10% share | pool reaching poorest 20% | unspent, returned |
|---|---|---|---|---|---|---|
| no zakat | 1.00 | 0.994 / 0.994 | 0.000 | 1.00 | - | - |
| profile:all_public | 0.66 | 0.364 / 0.333 | 0.084 | 0.52 | 0.20 | 0.00 |
| profile:poor_only | 0.57 | 0.189 / 0.097 | 0.144 | 0.47 | 0.52 | 0.00 |
| profile:eight_equal | 0.62 | 0.279 / 0.233 | 0.113 | 0.49 | 0.36 | 0.00 |
| profile:seven_equal | 0.62 | 0.279 / 0.233 | 0.113 | 0.49 | 0.38 | 0.00 |
| profile:needy_first | 0.59 | 0.223 / 0.152 | 0.132 | 0.48 | 0.46 | 0.00 |
| only:fuqara | 0.56 | 0.158 / 0.000 | 0.149 | 0.47 | 1.00 | 0.00 |
| only:masakin | 0.64 | 0.350 / 0.349 | 0.091 | 0.49 | 0.04 | 0.22 |
| only:amilin | 0.80 | 0.617 / 0.611 | 0.001 | 0.60 | 0.20 | 0.00 |
| only:muallafat | 0.65 | 0.357 / 0.333 | 0.074 | 0.48 | 0.28 | 0.00 |
| only:riqab | 0.61 | 0.272 / 0.231 | 0.117 | 0.47 | 0.93 | 0.00 |
| only:gharimin | 0.63 | 0.317 / 0.276 | 0.099 | 0.50 | 0.32 | 0.57 |
| only:sabilillah | 0.66 | 0.364 / 0.333 | 0.084 | 0.52 | 0.20 | 0.00 |
| only:ibn_sabil | 0.75 | 0.526 / 0.515 | 0.010 | 0.56 | 0.20 | 0.00 |

## B. Administrative overhead (weight of the amilin; the rest split between poor and needy)

| amilin weight | Gini | poverty before zakat | bottom-50% share | pool reaching poorest 20% |
|---|---|---|---|---|
| 0.0 | 0.57 | 0.189 | 0.144 | 0.52 |
| 0.125 | 0.59 | 0.225 | 0.132 | 0.43 |
| 0.25 | 0.61 | 0.270 | 0.117 | 0.34 |
| 0.4 | 0.64 | 0.333 | 0.094 | 0.28 |

## C. Zakat rate vs poverty headcount before zakat arrives (Gini in brackets)

| zakat/yr | all_public | poor_only | needy_first | eight_equal |
|---|---|---|---|---|
| 0.005 | 0.688 (Gini 0.88) | 0.655 (Gini 0.87) | 0.659 (Gini 0.87) | 0.667 (Gini 0.87) |
| 0.010 | 0.570 (Gini 0.80) | 0.510 (Gini 0.78) | 0.516 (Gini 0.78) | 0.530 (Gini 0.79) |
| 0.025 | 0.364 (Gini 0.66) | 0.189 (Gini 0.57) | 0.223 (Gini 0.59) | 0.279 (Gini 0.62) |
| 0.050 | 0.194 (Gini 0.52) | 0.039 (Gini 0.40) | 0.048 (Gini 0.42) | 0.133 (Gini 0.49) |
| 0.100 | 0.067 (Gini 0.38) | 0.012 (Gini 0.30) | 0.013 (Gini 0.30) | 0.029 (Gini 0.34) |

## D. Critical rate (tail exponent = 1) by allocation

| allocation | z_c |
|---|---|
| all_public | 0.00438 |
| poor_only | 0.00430 |
| eight_equal | 0.00451 |

## E. Nisab sensitivity (poverty line fixed at 0.2 x mean)

| setting | Gini | poverty before zakat | bottom-50% share |
|---|---|---|---|
| nisab=0.3,poor_only | 0.57 | 0.176 | 0.143 |
| nisab=0.3,eight_equal | 0.61 | 0.266 | 0.116 |
| nisab=0.5,poor_only | 0.57 | 0.189 | 0.144 |
| nisab=0.5,eight_equal | 0.62 | 0.279 | 0.113 |
| nisab=1.0,poor_only | 0.60 | 0.272 | 0.116 |
| nisab=1.0,eight_equal | 0.63 | 0.322 | 0.098 |

## F. Toggles: default (poor + needy) and each optional category switched on, equal shares

| enabled | Gini | poverty before zakat | bottom-50% share | pool reaching poorest 20% |
|---|---|---|---|---|
| default: fuqara+masakin | 0.57 | 0.189 | 0.144 | 0.52 |
| +amilin | 0.63 | 0.304 | 0.105 | 0.30 |
| +muallafat | 0.59 | 0.241 | 0.127 | 0.38 |
| +riqab | 0.58 | 0.220 | 0.135 | 0.59 |
| +gharimin | 0.57 | 0.195 | 0.141 | 0.57 |
| +sabilillah | 0.60 | 0.254 | 0.122 | 0.34 |
| +ibn_sabil | 0.62 | 0.284 | 0.112 | 0.32 |
