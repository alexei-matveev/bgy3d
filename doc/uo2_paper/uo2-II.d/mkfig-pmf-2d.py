#!/usr/bin/python3

import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

beta = 1.6889 # kcal^-1
T = 1 / beta

# PETSC files are big endian ...
# g3 = np.fromfile ("a1.bin", count=128**3)

def p(path):
    return "../../../test-qm/uranyl/0w/" + path

# Parsing text takes lots of time:
try:
    g3 = np.load (p("g0,g128.npy"))
except:
    g3 = np.loadtxt (p("g0,g128.txt"))
    np.save (p("g0,g128.npy"), g3)

g3.shape = (128, 128, 128)

# Less than a quater of a plane:
g2 = g3[64:64+32, 63, 64:64+32]

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.linspace(0, 5, 32)
Y = np.linspace(0, 5, 32)
X, Y = np.meshgrid(X, Y)

# Z truncated at 5 kcal here:
Z = - T * np.log (np.maximum (g2, 5.0e-4))

surf = ax.plot_surface (X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                        vmin = -3,
                        vmax = +3,
                        shade = True,
                        linewidth=0.25, antialiased=True)
cset = ax.contour(X, Y, Z, zdir='z', offset=-5, cmap=cm.coolwarm)

ax.set_zlim (-5.0, 5.0)
ax.set_xlabel ('x, A')
ax.set_ylabel ('y, A')
ax.set_zlabel ('W(x,y), kcal/mol')

ax.zaxis.set_major_locator(LinearLocator(3))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))

fig.colorbar (surf, shrink=1/2) #, aspect=10)

# Manually chosen:
ax.view_init (elev=50, azim=-75)

# plt.show()
plt.savefig (sys.argv[1], transparent=True, bbox_inches='tight')
