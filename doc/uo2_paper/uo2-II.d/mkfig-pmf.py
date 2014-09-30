#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function, division
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

def opt (x0, f):
    from scipy.optimize import minimize
    res = minimize (f, x0)
    return res.x

def one (path, ymax=8, rmax=3.0, linewidth=1.0):
    cols = loadtxt (path)
    r = cols[:, 0]
    rspline = spline (range (size (r)), r)
    rfine = map (rspline, arange (0, size (r), 1.0 / 8))
    gO = cols[:, 1]
    def pp (g, color):
        G = spline (r, g)
        W = lambda x: -T * log (G (x))

        plot (rfine, map (W, rfine),
              "-", label="W(r)", color=color, linewidth=linewidth)
        ylim (-2, 4)
        xlim (1, 7)
        axhline(0.0, color="black")
        ra = opt (2.5, lambda x: -G(x))
        rb = opt (4.5, lambda x: -G(x))
        rc = opt (3.0, lambda x: +G(x))
        print (path)
        print ("ra=", ra, "W(ra)=", W(ra))
        print ("rb=", rb, "W(ra)=", W(rb))
        print ("rc=", rc, "W(ra)=", W(rc))

    pp (gO, color="red")

figure (figsize=(8, 4))
rmax = 3.2
ymax = None

# MD data:
def p1 (path, color="black"):
    cols = loadtxt (path)
    r = cols[:, 0]
    g = cols[:, 1]
    G = spline (r, g)
    W = lambda x: -T * log (G(x))
    plot (r, map (W, r), "--", label="MD", color=color)

    ra = opt (2.5, lambda x: -G(x))
    rb = opt (4.5, lambda x: -G(x))
    rc = opt (3.0, lambda x: +G(x))
    print (path)
    print ("ra=", ra, "W(ra)=", W(ra))
    print ("rb=", rb, "W(ra)=", W(rb))
    print ("rc=", rc, "W(ra)=", W(rc))

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
