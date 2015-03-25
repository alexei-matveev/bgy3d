#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from numpy import loadtxt, shape, asarray
from pylab import *
import matplotlib.pyplot as plt

fontsize = 10
color = "black"
linewidth = 2

source = ["mm_acetic_acid_O-Hw.dat",
	  "gp_acetic_acid_O-Hw.dat",
	  "aq_acetic_acid_O-Hw.dat",
	  "mm_acetic_acid_H-Ow.dat",
	  "gp_acetic_acid_H-Ow.dat",
	  "aq_acetic_acid_H-Ow.dat",]

data = np.asarray([loadtxt(_) for _ in source])


def plot_one (d, title):

    # MM
    plt.plot(d[0, :, 0], d[0, :, 1], ":", color=color, lw=linewidth)
    # GP
    plt.plot(d[1, :, 0], d[1, :, 1], "--", color=color, lw=linewidth)
    # AQ
    plt.plot(d[2, :, 0], d[2, :, 1], "-", color=color, lw=linewidth)

    plt.xlim ((1, 8))

    plt.ylabel ("g(r)")

    plt.title (title, loc="left")


figure(figsize=(8, 8))

subplot (2, 1, 1)
plot_one (data[0:3, :, :], "O-HW")
plt.ylim ((0, 1.4))

subplot (2, 1, 2)
plot_one (data[3:6, :, :], "H-OW")
plt.ylim ((0, 1.2))
plt.xlabel (u"r, Å")


# show()
savefig (sys.argv[1], transparent=True, bbox_inches='tight')
