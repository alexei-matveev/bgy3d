#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot  as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline

kcal = 0.0433641195867 # eV

# Each of these has 4 columns (q, E(q), dG(q), G(q)):
soft = [loadtxt(f) for f in ("soft,interp.txt", "soft,optim.txt")]
hard = [loadtxt(f) for f in ("hard,interp.txt", "hard,optim.txt")]

def one ():
    def pp (two, label, color):
        interp, optim = two
        r_int = interp[:, 0]
        g_int = interp[:, 3] / kcal
        r_opt = optim[:, 0]
        g_opt = optim[:, 3] / kcal
        plot (r_int, g_int, "--", color=color)
        plot (r_opt, g_opt, "o", color=color, label=label)

    pp (soft, label="Soft", color="red")
    pp (hard, label= "Hard", color="blue")

figure (figsize=(8, 4))
one ()
xlabel ("q, A")
ylabel ("G(q), kcal/mol")
legend ()
# show()
savefig (sys.argv[1], transparent=True, bbox_inches='tight')
