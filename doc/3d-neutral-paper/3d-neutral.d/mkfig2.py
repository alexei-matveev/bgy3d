#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot as plt

fontsize=10
dmin, dmax = -0.05, 1.2
xdmax = dmax - 0.4
kw = {"linewidth": 1.5,
      "markeredgewidth": 1.5,
      "markersize": 8}

d = loadtxt("dipole.dat")

figure(figsize=(4, 2 * 4.25))
fig, (ax1, ax2) = subplots (1, 2, sharex=True, sharey=True)
ax1.set_aspect (1)
ax2.set_aspect (1)

# MM : d[:, 3]
# GP : d[:, 4]
# ax1.set_title('Gas phase', loc="left")
# also fit the MM vs. GP
fit = polyfit(d[:, 4], d[:, 3], 1)
fit_fn0 = poly1d(fit)
txt0 = "y = %.2fx + %.3f" % tuple (fit)
one = [dmin, xdmax]
ax1.plot (one, fit_fn0 (one), '-', color="0.00", **kw)
ax1.plot(d[:, 4], d[:, 3], 'o', color="white", label=txt0,  **kw)
print (polyfit(d[:, 4], d[:, 3], 1))

#xlim ((dmin, xdmax))
#ylim ((dmin, dmax))

# ax1.set_xlabel(u"Dipole (GP), eÅ")
# setp (ax1.get_xticklabels(), visible=False)
ax1.set_ylabel(u"Dipole (MM), eÅ")

# AQ, PCM vs. GP
# more space between subplots
# subplots_adjust(hspace=0.3)

# Line fit in the background:
fit = polyfit(d[:, 4], d[:, 0], 1)
fit_fn1 = poly1d(fit)
txt1 = "y = %.2fx + %.3f" % tuple (fit)

fit = polyfit(d[:, 4], d[:, 2], 1)
fit_fn2 = poly1d(fit)
txt2 = "y = %.2fx + %.3f" % tuple (fit)

# Line fit across the box:
two = [dmin, xdmax]
ax2.plot (two, fit_fn1 (two), '-', color="0.50",
          # label=txt1,
          **kw)
ax2.plot (two, fit_fn2 (two), '-', color="0.00",
          # label=txt2,
          **kw)

# PCM: d[:, 0]
# AQ: d[:, 2]
# ax2.set_title('Aqueous phase', loc="left")
ax2.plot(d[:, 4], d[:, 0], 's', color="0.75", label="PCM, " + txt1,
         markeredgecolor="0.50", **kw)
ax2.plot(d[:, 4], d[:, 2], 'o', color="1.00", label="SCF, " + txt2,
         **kw)

#xlim ((dmin, xdmax))
#ylim ((dmin, dmax))

ax1.set_xlabel(u"Dipole (Gas Phase), eÅ")
ax2.set_xlabel(u"Dipole (Gas Phase), eÅ")
ax2.set_ylabel(u"Dipole (Aqueous Phase), eÅ")
ax1.legend (ncol=1, numpoints=1, frameon=False,
	fontsize=fontsize, loc="lower right")
ax2.legend (ncol=1, numpoints=1, frameon=False,
	fontsize=fontsize, loc="lower right")

xmin, xmax = xlim()
ymin, ymax = ylim()
ax1.text (xmin + 0.1, ymax - 0.05, 'a', fontsize=15)
ax2.text (xmin + 0.1, ymax - 0.05, 'b', fontsize=15)

ax1.locator_params (tight=True, nbins=8)

# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
