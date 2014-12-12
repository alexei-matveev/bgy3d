#!/usr/bin/python

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot as plt

fontsize = 10

source = ["mm.dat", "fed.dat", "qm.dat", "scf-qm.dat"]

data = [loadtxt(_) for _ in source]


def sub(ax, i, title, do_legend=False):
    d = data[i]
    plt.title (title, loc="left", fontsize=fontsize)
    plt.plot(d[:, 0], d[:, 1], 'rs')
    plt.plot(d[:, 0], d[:, 2], 'bo')

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

    if i < 3:
	setp (ax.get_xticklabels(), visible=False)

figure(figsize=(4, 4 * 3))
sub(subplot(4, 1, 1), 0, "MM")
plt.ylabel(r'$\mathrm{\mathsf{\Delta G_{RISM} - \Delta G_{Expt}}}$')

sub(subplot(4, 1, 2), 1, "FED")
plt.ylabel(r'$\mathrm{\mathsf{\Delta G_{RISM} - \Delta G_{Expt}}}$')

sub(subplot(4, 1, 3), 2, "QM")
plt.ylabel(r'$\mathrm{\mathsf{\Delta G_{RISM} - \Delta G_{Expt}}}$')

sub(subplot(4, 1, 4), 3, "SCF QM")
plt.ylabel(r'$\mathrm{\mathsf{\Delta G_{rism} - \Delta G_{Expt}}}$')
plt.xlabel(r'$\mathrm{\mathsf{\rho * V}}$')

# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
