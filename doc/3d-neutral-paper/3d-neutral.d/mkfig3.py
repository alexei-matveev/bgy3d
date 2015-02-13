#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot as plt

fontsize = 13
linewidth = 1.5
emax = 4.1

d = loadtxt ("relax.dat")

figure (figsize=(4, 4))
ax = subplot(1, 1, 1)

# Place the line fit below the data points --- plot it first:
fit = polyfit (d[:, 0], d[:, 1], 1)
fn = poly1d (fit)

# Equation in the legend. Intercept happens to be negative:
a, b = fit
txt = u"δμ = %.2f δE - %.2f" % (a, -b)

# Plot line by two points:
plt.plot([0.0, emax], fn([0.0, emax]), '-',
         color="black", linewidth=linewidth,
         label = txt)

plt.plot(d[:, 0], d[:, 1], "o",
         color="white", markersize=8, markeredgewidth=linewidth)

#ax.set_title (u"δμ vs. δE", loc="left")
ax.set_xlabel (u"δE, kcal/mol")
ax.set_ylabel (u"δμ, kcal/mol")

xlim ((0.0, emax))
ylim ((-2 * emax, 0))
#plt.text(0.72, 0.75, 'a', fontsize=15)

ax.legend (ncol=1, numpoints=1, frameon=False,
           fontsize=fontsize, loc='lower left')
# plt.text(1.08, 1.12, 'b', fontsize=15)

ax.locator_params (tight=True, nbins=6)
# ax.xaxis.set_ticks (range(1, 4))

# show()
savefig (sys.argv[1], transparent=True, bbox_inches='tight')
