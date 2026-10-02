"""Values read off the figures of Kato (2022) by eye, for comparison only.

The paper publishes no tables and no code ("All relevant data are within the
paper"), so these are digitised from the printed plots and carry roughly
+/-0.02 reading error on top of Monte-Carlo noise. Treat `TOL` as the
tolerance for "reproduces the figure".
"""

TOL = 0.05

FIG4_TIMES = (10**3, 3 * 10**3, 10**4, 3 * 10**4, 10**5, 3 * 10**5, 10**6)

# Fig 4 / Fig 5 x-axis: Gini vs t. Keys are series labels.
# NB: the legend leader lines in Fig 4 cross: "J dw=0.2" labels the navy curve
# and "L rho=0" the pink one (checked against the 300 dpi render).
FIG4 = {
    "R":          (.27, .36, .39, .36, .36, .375, .38),
    "L rho=0.05": (.06, .11, .22, .48, .78, .92, .98),
    "L rho=0":    (.05, .09, .16, .27, .49, .75, .95),
    "J dw=0.2":   (.07, .12, .22, .37, .60, .86, .98),
    "J dw=0.1":   (.03, .06, .10, .19, .33, .54, .83),
}

# Fig 5 (xi = 0.5): only the curves that are legible. tp=1e4 curves coincide
# with exchange-only up to t=1e4 and then flatten.
FIG5 = {
    ("L", 10**5): (.06, .11, .22, .48, .78, .83, .84),
    ("J0.2", 10**5): (.07, .12, .22, .37, .58, .65, .66),
    ("J0.1", 10**5): (.03, .06, .11, .19, .34, .38, .38),
    ("R", 10**5): FIG4["R"],
}
FIG5_PLATEAU_TP1E4 = {"L": .31, "J0.2": .25, "J0.1": .13, "R": .38}

# Fig 6: Gini at the final time vs transfer rate xi.
FIG6_XI = (0.0, .125, .25, .375, .5, .75, 1.0)
FIG6 = {
    ("L", 10**5):    (.98, .91, .87, .85, .84, .81, .78),
    ("J0.2", 10**5): (.98, .80, .75, .69, .66, .62, .62),
    ("J0.1", 10**5): (.83, .58, .48, .41, .38, .36, .35),
    ("L", 10**4):    (.98, .64, .53, .40, .31, .24, .23),
    ("J0.2", 10**4): (.98, .42, .33, .27, .25, .23, .22),
    ("J0.1", 10**4): (.83, .22, .16, .14, .12, .12, .12),
}

# Fig 7 (xi=.5, tp=1e5): Gini vs q = 1..10 where the threshold is m_max / q.
FIG7_Q = tuple(range(1, 11))
FIG7 = {
    "L":    (.98, .89, .86, .84, .85, .84, .84, .84, .84, .84),
    "J0.2": (.98, .85, .80, .78, .75, .72, .73, .73, .75, .77),
    "J0.1": (.81, .65, .53, .41, .37, .37, .36, .37, .38, .47),
}
FIG7_T_LEVELS = {"L": .84, "J0.2": .66, "J0.1": .38}
