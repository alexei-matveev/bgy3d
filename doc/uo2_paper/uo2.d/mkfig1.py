#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot  as plt

fontsize = 14

# Each of these  has at least 1 +  3x3 columns for r and  all 3x3 pair
# distributions of uranyl and water sites:
src = [("cspce-pm13_rdf.dat", "cSPC/E", "-", "m"),
       ("shspce-pm13_rdf.dat", "SH-SPC/E", "--", "c"),
       ("pm13_rdf.dat", "PR-SPC/E", "-", "b")]

data, legends, styles, colors = zip(*src)
data = [loadtxt(_) for _ in data]

# Just these five columns, MD data does not have more than that:
cols = [0, 1, 2, 4, 5]
data = [d[:, cols] for d in data]

# Each of these has three columns (r, g(r), N(r)):
md = ["md-uo.dat", "md-uh.dat", "md-oo.dat", "md-oh.dat"]
md = [loadtxt(_) for _ in md]

# r-column is common for all of MD data, verify:
r = md[0][:, 0]
for d in md:
    assert (d[:, 0] == r).all()

# Five columns here, note the transposition:
md = array([r] + [d[:, 1] for d in md]).T

# To prepend MD entry zip them together again:
src = zip(data, legends, styles, colors)

src = [(md, "MD", ":", "k")] + src

data, legends, styles, colors = zip(*src)

figure()

def sub(ax, i, do_legend=False):

    for d, legend, style, color in zip(data, legends, styles, colors):
        plt.plot(d[:, 0], d[:, i],
                 linestyle=style, label=legend, lw=2, color=color,
                 dash_capstyle="round")

    # x0, x1 = plt.xlim()
    plt.xlim((1.0, 7.0))
    # y0, y1 = plt.ylim()
    plt.ylim(ymin=0.0)

    if do_legend:
        plt.legend(numpoints=1, frameon=False)

sub(subplot(221), 1, do_legend=True) # U-OW
plt.ylabel("g(r)", fontsize=fontsize)
plt.title("U-OW", fontsize=fontsize)

sub(subplot(222), 2)            # U-HW
plt.title("U-HW", fontsize=fontsize)

sub(subplot(223), 3)            # O-OW
plt.ylabel("g(r)", fontsize=fontsize)
plt.xlabel(u"r, Å", fontsize=fontsize)
plt.title("O-OW", fontsize=fontsize)

sub(subplot(224), 4)            # O-HW
plt.xlabel(u"r, Å", fontsize=fontsize)
plt.title("O-HW", fontsize=fontsize)

# show()
savefig(sys.argv[1], transparent=True)
