#!/usr/bin/python
# -*- coding: utf-8 -*-

from numpy import loadtxt
from pylab import *
import matplotlib.pyplot  as plt

fontsize = 14

source = [("gw96_rdf.dat", "GW", ":", "k"),
          ("kl1_rdf.dat", "KL1", "--", "m"),
          ("kl2_rdf.dat", "KL2", "--", "c"),
          ("spc_rdf.dat", "RM12", "-", "r"),
          ("pm13_rdf.dat", "PM13", "-", "b")]

data, legends, styles, colors = zip(*source)
data = [loadtxt(_) for _ in data]

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

sub(subplot(223), 4)            # O-OW
plt.ylabel("g(r)", fontsize=fontsize)
plt.xlabel(u"r, Å", fontsize=fontsize)
plt.title("O-OW", fontsize=fontsize)

sub(subplot(224), 5)            # O-HW
plt.xlabel(u"r, Å", fontsize=fontsize)
plt.title("O-HW", fontsize=fontsize)

# show()
savefig("test.svg", transparent=True)
