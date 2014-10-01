#!/usr/bin/python
from __future__ import division, print_function
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

N = 128
L = 10.0
H = 2 * L / N
# off-by-one on an 2L interval:
L1 = L * (N - 2) / N
def C (i):
    return i * H - L
assert L1 - C(N - 1) == 0.0

# Parsing text takes lots of time:
try:
    g3 = np.load (p("g0,g128.npy"))
except:
    g3 = np.loadtxt (p("g0,g128.txt"))
    np.save (p("g0,g128.npy"), g3)

g3.shape = (N, N, N)

ex = []
def extrema (g):
    # plt.imshow (g, extent=(-L, L, -L, L)) # cmap=cm.RdBu)
    # plt.show ()
    # x = range (N)
    # y = range (N)
    x = np.linspace (-L, L1, N)
    y = np.linspace (-L, L1, N)
    print (x)
    # x, y = np.meshgrid (x, y)
    from scipy.interpolate import RectBivariateSpline as spline
    from pts.func import NumDiff
    G = spline (x, y, g)
    print (g[0, 0], G(-L, -L))

    def f (r):
        x, y = r
        return G(x, y)

    f = NumDiff (f)
    def gn (r):
        df = f.fprime (r)
        return np.linalg.norm (df)

    from scipy.optimize import minimize

    def opt (f, x):
        res = minimize (f, x, options={"disp": True})
        return res.x

    xp = opt (lambda r: f(r), [0.0, 3.8])
    print ("x=", xp, "G=", f(xp), "|grad G|=", gn (xp))
    # ex.append (xp)
    xm = opt (lambda r: -f(r), [0.0, 2.4])
    print ("x=", xm, "G=", f(xm), "|grad G|=", gn (xm))
    # ex.append (xm)
    # xt = opt (lambda r: gn(r), [1.5, 2.95])
    # print ("x=", xt, "G=", f(xt), "|grad G|=", gn (xt))
    # ex.append (xt)
    xt = opt (lambda r: gn(r), [1.5, 2.8])
    print ("x=", xt, "G=", f(xt), "|grad G|=", gn (xt))
    ex.append (xt)

n = N // 2
extrema (g3[:, n, :])
ex = np.array (ex)
print (ex)

# Less than a quater of a plane:
g2 = g3[n:n + n // 2, n, n:n + n // 2]

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.linspace(C(n), C(n + n // 2 - 1), n // 2)
Y = np.linspace(C(n), C(n + n // 2 - 1), n // 2)
X, Y = np.meshgrid(X, Y)

# Z truncated at 5 kcal here:
Z = - T * np.log (np.maximum (g2, 5.0e-4))

surf = ax.plot_surface (X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                        vmin = -3,
                        vmax = +3,
                        shade = True,
                        linewidth=0.25, antialiased=True)
cset = ax.contour(X, Y, Z, zdir='z', offset=-5, cmap=cm.coolwarm)
ax.scatter (ex[:, 1], ex[:, 0], zs=-5, c="black")

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
