#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function

import sys
from numpy import loadtxt, array, shape, vstack, copy, min, max
from pylab import *
import matplotlib.pyplot  as plt
from scipy.interpolate import UnivariateSpline as spline

colorS = "#c0c0c0"              # silver
colorH = "#c0c0c0"
colorQ = "black" # "#000000"              # black

kcal = 0.0433641195867 # eV

# Each of these has 4 columns (q, E(q), dG(q), G(q)):
soft = [loadtxt(f) for f in ("soft,interp.txt", "soft,optim.txt")]
hard = [loadtxt(f) for f in ("hard,interp.txt", "hard,optim.txt")]
qmrs, e0 = loadtxt ("qmrism,stationary.txt"), -777424.35579213456409218810

if False:
    d, e0 = loadtxt ("qmrism,optim.txt"), -777423.62550623293672912596
    print ("shape(d)=", shape (d))
    d0 = copy (d[[0, 0]])
    d0[0, 0] = - d0[0, 0]
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
    def pp (two, symbol, label, color):
        interp, optim = two
        r_int = interp[:, 0]
        g_int = interp[:, 3] / kcal
        r_opt = optim[:, 0]
        g_opt = optim[:, 3] / kcal
        plot (r_int, g_int, "--", color=color)
        plot (r_opt, g_opt, symbol, color=color, label=label)

    # For markers see http://matplotlib.org/api/markers_api.html
    pp (soft, "s", label="Soft", color=colorS)
    pp (hard, "d", label= "Hard", color=colorH)
    def pq ():
        r = qmrs[:, 0]
        e = (qmrs[:, 3] - e0) / kcal
        print ("r=", r)
        plot (r, e, "o-", label="QM", color=colorQ)

        # espline = spline (r, e, s=0.0)
        # rfine = arange (min (r), max (r), 0.02)
        # plot (rfine, map (espline, rfine), "--", color=colorQ)

    pq()


figure (figsize=(6, 6))
one ()
xlabel ("q, A")
xlim ((-0.02, 2.0))
ylabel ("G(q), kcal/mol")
ylim ((-435., -423.))
legend (loc=3, numpoints=1)
show()
#savefig (sys.argv[1], transparent=True, bbox_inches='tight')
