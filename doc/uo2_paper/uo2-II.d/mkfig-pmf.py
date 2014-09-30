#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape, log
from pylab import *
import matplotlib.pyplot  as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline

def p (path):
    return "/home/matveev/darcs/bgy3d/test-qm/uranyl/sym/rdfs-10x256-trilin/" + path

# Each of these has three columns (r, g(r), N(r)):
md = ["md-uo.dat", "md-uh.dat", "md-oo.dat", "md-oh.dat"]
md = [loadtxt("../uo2.d/" + _) for _ in md]

rho = 0.0333295 # A^-3
beta = 1.6889 # kcal^-1
T = 1 / beta

def one (path, ymax=8, rmax=3.0, linewidth=1.0):
    cols = loadtxt (path)
    r = cols[:, 0]
    rspline = spline (range (size (r)), r)
    rfine = map (rspline, arange (0, size (r), 1.0 / 8))
    gO = cols[:, 1]
    def pp (g, color):
        gspline = spline (r, g)

        plot (rfine , -T * log (map (gspline, rfine)),
              "-", label="g(r)", color=color, linewidth=linewidth)
        ylim (-2, 4)
        xlim (1, 7)
        axhline(0.0, color="black")

    pp (gO, color="red")

figure (figsize=(8, 4))
rmax = 3.2
ymax = None

# MD data:
def p1 (path, color="black"):
    cols = loadtxt (path)
    r = cols[:, 0]
    g = cols[:, 1]
    plot (r, -T * log (g), "--", label="MD", color=color)

def p2 (xo, xh, x1, x2, x3):
    # FIXME: MD data uses PM13 uranyl:
    if True:
        p1 (xo, color="black")

    # RISM data:
    one (p (x1), ymax=ymax, rmax=rmax, linewidth=1.0)
    one (p (x2), ymax=ymax, rmax=rmax, linewidth=1.6)
    one (p (x3), ymax=ymax, rmax=rmax, linewidth=2.56)

p2 ("../uo2.d/md-uo.dat",
    "../uo2.d/md-uh.dat",
    "0w,PSE3,U.rdf",
    "0w,PSE2,U.rdf",
    "0w,KH,U.rdf")
# title ("PMF")
ylabel ("W(r), kcal/mol")
xlabel ("r, A")

# show()
savefig (sys.argv[1], transparent=True, bbox_inches='tight')
