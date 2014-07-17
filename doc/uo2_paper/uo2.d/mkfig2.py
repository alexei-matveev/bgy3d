#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt
from pylab import *
import matplotlib.pyplot  as plt

fontsize = 10

source = [("gw96_rdf.dat", "GW", ":", "k"),
          ("kl1_rdf.dat", "KL1", "--", "m"),
          ("kl2_rdf.dat", "KL2", "--", "c"),
          ("spc_rdf.dat", "RM12", "-", "r"),
          ("pm13_rdf.dat", "PM13", "-", "b")]

data, legends, styles, colors = zip(*source)
data = [loadtxt(_) for _ in data]

# Better be longer for a 4 x 1 subplot grid. Dimensions in inches:
figure(figsize=(4, 4 * 1.5))

def sub(ax, i, title, do_legend=False):

    plt.title (title, loc="left", fontsize=fontsize)
    for d, legend, style, color in zip(data, legends, styles, colors):
        plt.plot(d[:, 0], d[:, i],
                 linestyle=style, label=legend, lw=2, color=color,
                 dash_capstyle="round")

    # x0, x1 = plt.xlim()
    plt.xlim((1.0, 7.0))
    # y0, y1 = plt.ylim()
    plt.ylim(ymin=0.0)
    ax.locator_params (tight=True, nbins=6)

    if do_legend:
        plt.legend (ncol=2, numpoints=1, frameon=False,
                    fontsize=fontsize)

    if i != 5:
        setp (ax.get_xticklabels(), visible=False)

sub(subplot(4, 1, 1), 1, "U-OW", do_legend=True)
plt.ylabel("g(r)", fontsize=fontsize)

sub(subplot(4, 1, 2), 2, "U-HW")
plt.ylabel("g(r)", fontsize=fontsize)

sub(subplot(4, 1, 3), 4, "O-OW")
plt.ylabel("g(r)", fontsize=fontsize)

sub(subplot(4, 1, 4), 5, "O-HW")
plt.ylabel("g(r)", fontsize=fontsize)
plt.xlabel(u"r, Å", fontsize=fontsize)

# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
