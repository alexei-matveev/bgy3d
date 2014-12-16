#!/usr/bin/python

import sys
from numpy import loadtxt, array, shape
from pylab import *
import matplotlib.pyplot as plt

fontsize=10

d = loadtxt("dipole.dat")

figure(figsize=(4, 2 * 4.25))

# MM and GP
ax1 = subplot(2, 1, 1)
ax1.set_title('Gas phase', loc="left")
plt.plot(d[:, 3], d[:, 4], 'rs')
ax1.set_xlabel('MM')
ax1.set_ylabel('FED')

# QM, SCF-QM vs. PCM
ax2 = subplot(2, 1, 2)
ax2.set_title('Aqueous phase', loc="left")
plt.plot(d[:, 0], d[:, 1], 'rs', label='QM')
plt.plot(d[:, 0], d[:, 2], 'bo', label='SCF QM')

# line fit
linewidth=2.0
fit = polyfit(d[:, 0], d[:, 1], 1)
fit_fn1 = poly1d(fit)
fit = polyfit(d[:, 0], d[:, 2], 1)
fit_fn2 = poly1d(fit)

plt.plot(d[:, 0], fit_fn1(d[:, 0]), 'r--', linewidth=linewidth)
plt.plot(d[:, 0], fit_fn2(d[:, 0]), 'b--', linewidth=linewidth)

ax2.set_xlabel('PCM')
ax2.set_ylabel('RISM')
ax2.legend (ncol=1, numpoints=1, frameon=False,
	fontsize=fontsize, loc='upper left')

# show()
savefig(sys.argv[1], transparent=True, bbox_inches='tight')
