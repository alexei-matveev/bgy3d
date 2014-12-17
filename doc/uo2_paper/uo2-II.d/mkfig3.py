#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
from numpy import loadtxt, array, shape, vstack, copy, min, max
from pylab import *
import matplotlib.pyplot  as plt
from scipy.interpolate import UnivariateSpline as spline

colorH = "#000000"              # black
colorS = "#c0c0c0"              # silver

kcal = 0.0433641195867 # eV

# Each of these has 4 columns (q, E(q), dG(q), G(q)):
soft = [loadtxt(f) for f in ("soft,interp.txt", "soft,optim.txt")]
hard = [loadtxt(f) for f in ("hard,interp.txt", "hard,optim.txt")]
d, e0 = loadtxt ("qmrism,optim.txt"), -777423.62550623293672912596
print ("shape(d)=", shape (d))
d1 = copy (d[1:-2])             # (1.7 -> 0.2)
d2 = copy (d[1:-1])             # (1.7 -> 0)
print ("d1=", d1)
print ("shape(d1)=", shape (d1))
d1[:, 0] = - d1[:, 0]           # (-1.7 -> -0.2)
d2 = d2[::-1]                   # (0 -> 1.7)
print ("d1=", d1)
qmrs = vstack ([d1, d2])        # (-1.7 -> -0.2, 0 -> 1.7)
print ("shape(qmrs)=", shape (qmrs))
print ("d=", qmrs)

def one ():
    def pp (two, label, color):
        interp, optim = two
        r_int = interp[:, 0]
        g_int = interp[:, 3] / kcal
        r_opt = optim[:, 0]
        g_opt = optim[:, 3] / kcal
        plot (r_int, g_int, "--", color=color)
        plot (r_opt, g_opt, "o", color=color, label=label)

    pp (soft, label="Soft", color=colorS)
    pp (hard, label= "Hard", color=colorH)
    def pq ():
        r = qmrs[:, 0]
        e = (qmrs[:, 3] - e0) / kcal
        print ("r=", r)
        plot (r, e, "o", color="red")

        espline = spline (r, e, s=0.0)
        rfine = arange (min (r), max (r), 0.02)
        plot (rfine, map (espline, rfine), "--", color="red")

    pq()


figure (figsize=(8, 4))
one ()
xlabel ("q, A")
ylabel ("G(q), kcal/mol")
legend ()
# show()
savefig (sys.argv[1], transparent=True, bbox_inches='tight')
