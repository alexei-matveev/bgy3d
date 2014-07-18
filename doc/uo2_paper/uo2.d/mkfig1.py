#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot  as plt

fontsize = 10

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

# Better be longer for a 4 x 1 subplot grid. Dimensions in inches:
figure(figsize=(4, 4 * 2))

def sub(ax, i, title, do_legend=False):

    # ax.set(aspect="equal")
    plt.title (title, loc="left", fontsize=fontsize)
    for d, legend, style, color in zip(data, legends, styles, colors):
        plt.plot(d[:, 0], d[:, i],
                 linestyle=style, label=legend, lw=2, color=color,
                 dash_capstyle="round")

    # x0, x1 = plt.xlim()
    plt.xlim((1.0, 7.0))
    # y0, y1 = plt.ylim()
    plt.ylim(ymin=0.0)
    ax.locator_params (tight=True, nbins=3)

    if do_legend:
        plt.legend (ncol=1, numpoints=1, frameon=False,
                    fontsize=fontsize)

    if i < 4:
        setp (ax.get_xticklabels(), visible=False)

sub (subplot(4, 1, 1), 1, "U-OW", do_legend=True)
plt.ylabel("g(r)", fontsize=fontsize)

sub (subplot(4, 1, 2), 2, "U-HW")
plt.ylabel("g(r)", fontsize=fontsize)

sub (subplot(4, 1, 3), 3, "O-OW")
plt.ylabel("g(r)", fontsize=fontsize)

sub (subplot(4, 1, 4), 4, "O-HW")
plt.ylabel("g(r)", fontsize=fontsize)
plt.xlabel(u"r, Å", fontsize=fontsize)

# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
