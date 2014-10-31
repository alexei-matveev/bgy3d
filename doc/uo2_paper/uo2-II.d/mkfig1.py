#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot  as plt
from scipy.interpolate import InterpolatedUnivariateSpline as spline

colorOW = "#000000"              # black
colorHW = "#C0C0C0"              # silver

def p (path):
    return "/home/matveev/darcs/bgy3d/test-qm/uranyl/sym/rdfs-10x256-trilin/" + path

# Each of these has three columns (r, g(r), N(r)):
md = ["md-uo.dat", "md-uh.dat", "md-oo.dat", "md-oh.dat"]
md = [loadtxt("../uo2.d/" + _) for _ in md]

rho = 0.0333295 # A^-3

def one (path, ymax=8, rmax=3.0, linewidth=1.0):
    cols = loadtxt (path)
    r = cols[:, 0]
    rspline = spline (range (size (r)), r)
    rfine = map (rspline, arange (0, size (r), 1.0 / 8))
    gO = cols[:, 1]
    gH = cols[:, 2]
    def pp (g, color):
        gspline = spline (r, g)
        Nprime = spline (r, 4 * pi * rho * r**2 * g)
        N = Nprime.antiderivative()

        plot (rfine, map (gspline, rfine), "-", color="white", linewidth=(linewidth + 2.0))
        plot (rfine , map (gspline, rfine), "-", label="g(r)", color=color, linewidth=linewidth)
        # plot (r, map (N, r), "--", label="n(r)", color=color)
        # legend ()
        ylim (0) #, ymax)
        xlim (1, 7)
        # xlabel ("r, A")
        # ylabel ("")
        # title (path)
        # axvline(rmax, color="black")
        # axhline(1.0, color="black")
        print "n(", rmax, ") =", N(rmax), color, path
        # print "n(", r[-1], ") =", N(r[-1])

    pp (gO, color=colorOW)
    pp (gH, color=colorHW)

figure (figsize=(8, 6))
rmax = 3.2
ymax = None

# MD data:
def p1 (path, color="black"):
    cols = loadtxt (path)
    r = cols[:, 0]
    g = cols[:, 1]
    Nprime = spline (r, 4 * pi * rho * r**2 * g)
    N = Nprime.antiderivative()
    print "n(", rmax, ") =", N(rmax), "(MD)"
    plot (r, g, "--", label="g(r)", color=color)

def p2 (xo, xh, x1, x2, x3):
    # FIXME: MD data uses PM13 uranyl:
    if False:
        p1 (xo, color=colorOW)
        p1 (xh, color=colorHW)

    # RISM data:
    one (p (x1), ymax=ymax, rmax=rmax, linewidth=1.0)
    one (p (x2), ymax=ymax, rmax=rmax, linewidth=1.6)
    one (p (x3), ymax=ymax, rmax=rmax, linewidth=2.56)

#
# For another figure replace 0 <-> 5:
#
subplot (2, 1, 1)
p2 ("../uo2.d/md-uo.dat",
    "../uo2.d/md-uh.dat",
    "5w,PSE3,U.rdf",
    "5w,PSE2,U.rdf",
    "5w,KH,U.rdf")
title ("U-OW, U-HW", loc="left")
ylabel ("g(r)")

subplot (2, 1, 2)
p2 ("../uo2.d/md-oo.dat",
    "../uo2.d/md-oh.dat",
    "5w,PSE3,OU.rdf",
    "5w,PSE2,OU.rdf",
    "5w,KH,OU.rdf")
title ("O-OW, O-HW", loc="left")
ylabel ("g(r)")
xlabel ("r, A")

# show()
savefig (sys.argv[1], transparent=True, bbox_inches='tight')
