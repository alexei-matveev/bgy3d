#!/usr/bin/python

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot as plt

fontsize=10
dmin, dmax = -0.05, 1.2
xpos, ypos = dmax - 0.12, dmin + 0.08
kw = {"linewidth": 1.5,
      "markeredgewidth": 1.5,
      "markersize": 8}

d = loadtxt("dipole.dat")

figure(figsize=(4, 2 * 4.25))

# MM and GP
ax1 = subplot(2, 1, 1)
# ax1.set_title('Gas phase', loc="left")
plt.plot(d[:, 3], d[:, 4], 'o', color="white", **kw)

xlim ((dmin, dmax))
ylim ((dmin, dmax))

ax1.set_xlabel('Dipole (MM), Debye')
ax1.set_ylabel('Dipole (FED), Debye')
plt.text (xpos, ypos, 'a', fontsize=15)

# QM, SCF-QM vs. PCM
ax2 = subplot(2, 1, 2)

# more space between subplots
subplots_adjust(hspace=0.3)

# Line fit in the background:
fit = polyfit(d[:, 0], d[:, 1], 1)
fit_fn1 = poly1d(fit)
txt1 = "y = %.2fx + %.2f" % tuple (fit)

fit = polyfit(d[:, 0], d[:, 2], 1)
fit_fn2 = poly1d(fit)
txt2 = "y = %.2fx + %.2f" % tuple (fit)

# Line fit across the box:
two = [dmin, dmax]
plt.plot (two, fit_fn1 (two), '-', color="0.75",
          # label=txt1,
          **kw)
plt.plot (two, fit_fn2 (two), '-', color="0.00",
          # label=txt2,
          **kw)

# ax2.set_title('Aqueous phase', loc="left")
plt.plot(d[:, 0], d[:, 1], 's', color="0.75", label="QM, " + txt1,
         markeredgecolor="0.50", **kw)
plt.plot(d[:, 0], d[:, 2], 'o', color="1.00", label="SCF-QM, " + txt2,
         **kw)

xlim ((dmin, dmax))
ylim ((dmin, dmax))

ax2.set_xlabel('Dipole (PCM), Debye')
ax2.set_ylabel('Dipole (RISM), Debye')
ax2.legend (ncol=1, numpoints=1, frameon=False,
	fontsize=fontsize, loc="upper left")
plt.text (xpos, ypos, 'b', fontsize=15)


# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
