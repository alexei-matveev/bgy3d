#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape, arange
from pylab import *
import matplotlib.pyplot as plt

fontsize = 10
kw = {"linewidth": 1.5,
      "markeredgewidth": 1.0,
      "markersize": 6}


source = ["mm.dat", "pt1.dat", "pt2.dat", "scf.dat"]

data = [loadtxt(_) for _ in source]


def sub(ax, i, title, do_legend=False):
    plt.axhline (color="black")
    d = data[i]

    # Linear fit in the background:
    fit = polyfit(d[:, 0], d[:, 1], 1)
    fit_fn1 = poly1d(fit)
    fit = polyfit(d[:, 0], d[:, 2], 1)
    fit_fn2 = poly1d(fit)

    two = [1.6, 5.6]
    plt.plot (two, fit_fn1 (two), '-', color="0.50", **kw)
    plt.plot (two, fit_fn2 (two), '-', color="0.00", **kw)

    # plt.title (title, loc="center", fontsize=fontsize)
    plt.plot(d[:, 0], d[:, 1], 's', label='uncorrected', color="0.50", mfc="0.75", mec="0.50", **kw)
    plt.plot(d[:, 0], d[:, 2], 'o', label='corrected', color="1.00", mfc="1.00", **kw)

    # xlim
    plt.xlim((1.5, 6.0))

    plt.ylim((-10.0, 30.0))

    # ticks
    xticks = np.arange (2.0, 6.1, 1.0)
    yticks = np.arange (-10, 31, 10)
    plt.xticks (xticks)
    plt.yticks (yticks)

    # text label in plot
    plt.text (5.7, 26, title, fontsize = 15)

    if do_legend:
	plt.legend (ncol=1, numpoints=1, frameon=False,
		fontsize=fontsize, loc='upper left')

    if i < 3:
	setp (ax.get_xticklabels(), visible=False)

figure(figsize=(4, 4 * 3))
# MM
sub(subplot(4, 1, 1), 0, "a", do_legend=False)
# plt.ylabel(r'$\mathrm{\mathsf{\Delta^2 G}}$', rotation='horizontal', )
plt.text (1.2, 33, u"Δ$^\mathregular{2}$G, kcal/mol", rotation='horizontal')

# FED
sub(subplot(4, 1, 2), 1, "b")

# QM
sub(subplot(4, 1, 3), 2, "c")

# SCF QM
sub(subplot(4, 1, 4), 3, "d")
plt.xlabel (u"ρV")

# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
