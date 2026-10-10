"""Values read off the figures of Kato (2022) for comparison only.

The paper publishes no tables and no code ("All relevant data are within the
paper"). Series marked [digitised] were extracted programmatically from the 305 dpi
raster images embedded in the PDF (experiments/digitize_kato_figures.py, y-axis
calibrated on the tick labels, ~0.002 Gini per pixel); they replaced earlier by-eye
readings that agreed to within 0.010 (mean 0.004). Series without the tag are still by
eye (markers hidden behind other markers, or hollow markers) at roughly +/-0.02. The
published points are single runs; Kato's own run-to-run scatter is visible in Fig 6, where
the xi=0 point of one series sits ~0.03 away from the same configuration in another
series. Treat `TOL` as the tolerance for "reproduces the figure".
"""

TOL = 0.05

FIG4_TIMES = (10**3, 3 * 10**3, 10**4, 3 * 10**4, 10**5, 3 * 10**5, 10**6)

# Fig 4 / Fig 5 x-axis: Gini vs t. Keys are series labels.
# NB: the legend leader lines in Fig 4 cross: "J dw=0.2" labels the navy curve
# and "L rho=0" the pink one (checked against the 300 dpi render).
FIG4 = {
    "R":          (.27, .36, .39, .36, .36, .375, .38),
    "L rho=0.05": (.06, .11, .22, .480, .781, .922, .98),       # [digitised] t=3e4..3e5
    "L rho=0":    (.05, .089, .156, .276, .493, .756, .951),     # [digitised] except t=1e3
    "J dw=0.2":   (.067, .117, .218, .376, .597, .855, .988),   # [digitised]
    "J dw=0.1":   (.036, .060, .101, .189, .331, .543, .839),   # [digitised]
}

# Fig 5 (xi = 0.5): only the curves that are legible. tp=1e4 curves coincide
# with exchange-only up to t=1e4 and then flatten.
FIG5 = {
    ("L", 10**5): (.06, .11, .22, .483, .772, .838, .842),       # [digitised] t>=3e4
    ("J0.2", 10**5): (.065, .120, .224, .363, .578, .652, .661),  # [digitised]
    ("J0.1", 10**5): (.038, .063, .109, .191, .343, .380, .380),  # [digitised]
    ("R", 10**5): FIG4["R"],
}
FIG5_PLATEAU_TP1E4 = {"L": .31, "J0.2": .25, "J0.1": .13, "R": .38}

# Fig 6: Gini at the final time vs transfer rate xi.
FIG6_XI = (0.0, .125, .25, .375, .5, .75, 1.0)
FIG6 = {
    # tp=1e5, xi=0 point is the exchange-only value; xi>0 [digitised]
    ("L", 10**5):    (.98, .907, .870, .850, .842, .814, .782),
    ("J0.2", 10**5): (.98, .802, .749, .689, .661, .621, .621),
    ("J0.1", 10**5): (.83, .584, .480, .411, .379, .356, .343),
    ("L", 10**4):    (.98, .64, .53, .40, .31, .24, .23),
    ("J0.2", 10**4): (.98, .42, .33, .27, .25, .23, .22),
    ("J0.1", 10**4): (.83, .22, .16, .14, .12, .12, .12),
}

# Fig 7 (xi=.5, tp=1e5): Gini vs q = 1..10 where the threshold is m_max / q.
FIG7_Q = tuple(range(1, 11))
FIG7 = {
    # [digitised] (L at q=1 is hidden behind the J0.2 marker)
    "L":    (.98, .892, .859, .843, .851, .843, .847, .839, .839, .839),
    "J0.2": (.989, .855, .806, .782, .753, .721, .729, .733, .753, .774),
    "J0.1": (.810, .651, .533, .410, .370, .374, .367, .378, .390, .471),
}
FIG7_T_LEVELS = {"L": .84, "J0.2": .66, "J0.1": .38}

# Fig 3 (single run, t = 1e6): number of the 1000 agents in the first 0.1-wide wealth bin
# [0, 0.1), read from the native-resolution bar heights (+/-10). Wealth is in raw units, not
# normalised, so it also constrains total wealth: Kato's t=1e3/1e4 panels centre on m = 1.
FIG3_FIRST_BIN_1E6 = {"J dw=0.1": 390, "J dw=0.2": 910, "L rho=0.05": 885}
FIG3_TOL_AGENTS = 50
