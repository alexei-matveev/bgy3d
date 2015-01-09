#!/usr/bin/python

import sys
from numpy import loadtxt, array, shape, arange
from pylab import *
import matplotlib.pyplot as plt

fontsize = 10

source = ["mm.dat", "fed.dat", "qm.dat", "scf-qm.dat"]

data = [loadtxt(_) for _ in source]


def sub(ax, i, title, do_legend=False):
    d = data[i]
    # plt.title (title, loc="center", fontsize=fontsize)
    plt.plot(d[:, 0], d[:, 1], 'rs', label='uncorrected')
    plt.plot(d[:, 0], d[:, 2], 'bo', label='corrected')

    # linear fit
    linewidth=2.0
    fit = polyfit(d[:, 0], d[:, 1], 1)
    fit_fn1 = poly1d(fit)
    fit = polyfit(d[:, 0], d[:, 2], 1)
    fit_fn2 = poly1d(fit) 


    plt.plot(d[:, 0], fit_fn1(d[:, 0]), 'r--', linewidth=linewidth)
    plt.plot(d[:, 0], fit_fn2(d[:, 0]), 'b--', linewidth=linewidth)

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
plt.text (1.2, 33, r'$\mathrm{\mathsf{\Delta^{2}G}}$', rotation='horizontal')

# FED
sub(subplot(4, 1, 2), 1, "b")

# QM
sub(subplot(4, 1, 3), 2, "c")

# SCF QM
sub(subplot(4, 1, 4), 3, "d")
plt.xlabel(r'$\mathrm{\mathsf{\rho V}}$')

# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
