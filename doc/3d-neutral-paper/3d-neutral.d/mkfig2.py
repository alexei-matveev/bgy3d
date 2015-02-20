#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot as plt

fontsize=10
dmin, dmax = -0.05, 1.2
xdmax = dmax - 0.4
xpos, ypos = xdmax - 0.12, dmin + 0.08
kw = {"linewidth": 1.5,
      "markeredgewidth": 1.5,
      "markersize": 8}

d = loadtxt("dipole.dat")

figure(figsize=(4, 2 * 4.25))

# MM : d[:, 3]
# GP : d[:, 4]
ax1 = subplot(2, 1, 1)
# ax1.set_title('Gas phase', loc="left")
plt.plot(d[:, 4], d[:, 3], 'o', color="white", **kw)

xlim ((dmin, xdmax))
ylim ((dmin, dmax))

# ax1.set_xlabel(u"Dipole (GP), eÅ")
setp (ax1.get_xticklabels(), visible=False)
ax1.set_ylabel(u"Dipole (MM), eÅ")
plt.text (xpos, ypos, 'a', fontsize=15)

# AQ, PCM vs. GP
ax2 = subplot(2, 1, 2)

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
plt.plot (two, fit_fn1 (two), '-', color="0.75",
          # label=txt1,
          **kw)
plt.plot (two, fit_fn2 (two), '-', color="0.00",
          # label=txt2,
          **kw)

# PCM: d[:, 0]
# AQ: d[:, 2]
# ax2.set_title('Aqueous phase', loc="left")
plt.plot(d[:, 4], d[:, 0], 's', color="0.75", label="PCM, " + txt1,
         markeredgecolor="0.50", **kw)
plt.plot(d[:, 4], d[:, 2], 'o', color="1.00", label="AQ, " + txt2,
         **kw)

xlim ((dmin, xdmax))
ylim ((dmin, dmax))

ax2.set_xlabel(u"Dipole (Gas Phase), eÅ")
ax2.set_ylabel(u"Dipole (Aqueous Phase), eÅ")
ax2.legend (ncol=1, numpoints=1, frameon=False,
	fontsize=fontsize, loc="upper left")
plt.text (xpos, ypos, 'b', fontsize=15)


# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
