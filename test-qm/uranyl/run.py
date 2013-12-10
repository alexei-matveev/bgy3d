from __future__ import with_statement
#
# Tell Python interpreter where to find custom modules:
#
#   export PYTHONPATH=~/PYTHON:~/darcs/bgy3d/python
#
# For the  rubgy.scm script to run  one might need to  set the library
# search path:
#
# export LD_LIBRARY_PATH=~/darcs/bgy3d
#
import os
from ase.calculators.paragauss import ParaGauss
from ase.io import read
from pts.cfunc import Affine, Cartesian
from pts.memoize import Memoize, DirStore
from pts.units import kcal, eV, Hartree
from pts.func import compose, NumDiff
from pts.zmat import ZMat
from pts.qfunc import QFunc
from pts.fopt import minimize
from rism import Server

#
# FIXME: PG will  choke on gxfile, and saved_scfstate  from a previous
# calculation unless you clean:
#
def clean ():
    files = ["./gxfile", "./saved_scfstate.dat"]
    for path in files:
        try:
            os.unlink (path)
        except:
            pass

ucmd = \
"""
/users/alexei/darcs/bgy3d-wheezy/guile/runbgy.scm
--solvent "water, SPC/E" --rho 0.0333295 --beta 1.6889 --L 160 --N 4096
--solute "uranyl, SPC" --norm-tol 1e-14 --dielectric 78.4
"""

wcmd = \
"""
/users/alexei/darcs/bgy3d-wheezy/guile/runbgy.scm
--solvent "water, SPC/E" --rho 0.0333295 --beta 1.6889 --L 160 --N 4096
--solute "water, SPC/E" --norm-tol 1e-14 --dielectric 78.4
"""

command = "salloc -n8 mpirun ~/darcs/ttfs-mac/runqm"

ucalc = ParaGauss (cmdline=command, input="uranyl.scm")
wcalc = ParaGauss (cmdline=command, input="water.scm")

uatoms = read ("uo22+,gp.xyz")
watoms = read ("h2o.xyz")

x = uatoms.get_positions ()

if False:
    trafo = compose (Cartesian (), Affine (x.reshape (x.size, 1)))
else:
    # Z-matrix for a linear UOO molecule:
    zm = ZMat ([(None, None, None),
                (0, None, None),
                (0, 1, None)])

    # 2 -> 3 linear trafo
    reduc = Affine ([[1, 0],
                     [1, 0],
                     [0, 1]])

    # Two bond length are constrained to be the same:
    trafo = compose (zm, reduc)

# Internal coordinates:
s = trafo.pinv (x)
print "XXX: internals = ", s
print  trafo (s)
print "XXX: fprime = ", trafo.fprime (s)
# exit (1)

clean ()

with QFunc (uatoms, ucalc) as f:
    with Server (ucmd) as g:
        f = Memoize (f, DirStore (salt="uo22+, gp, qm Dec 9"))
        g = Memoize (g, DirStore (salt=ucmd + "Dec 9"))

        # One could compose  (f + g, trafo), but we  use f(s) and g(s)
        # below as functions of internal coordinates:
        f = compose (f, trafo)
        g = compose (g, trafo)

        e = f + g

        print "s(start)=", s
        print "x(start)=\n", trafo (s)

        s, info = minimize (e, s, algo=1)
        print s, e (s), "eV", e (s) / kcal, "kcal", info["converged"]

        print "converged =", info["converged"], "in", info["iterations"], "iterations"
        print "s(min)=", s
        print "x(min)=\n", trafo (s)

        units = [(kcal, "kcal"), (eV, "eV"), (Hartree, "Hartree")]
        for u, uu in units:
            print "e = ", f(s)/u, "+", g(s)/u, "=", e(s)/u, "(%s)" % uu
        for u, uu in units:
            print "g = ", f.fprime(s)/u, "+", g.fprime(s)/u, "=", e.fprime(s)/u, "(%s/Unit)" % uu



print "==================== WATER ======================"
x = watoms.get_positions ()

#
# Z-matrix for water:
#
zmt = ZMat ([(None, None, None),
             (1, None, None),
             (1, 2, None)], base=1)
trafo = zmt

#
# Initial values of internal coordinates:
#
s = trafo.pinv (x)
assert max (abs (s - trafo.pinv (zmt (s)))) < 1.0e-10

clean ()

with QFunc (watoms, wcalc) as f, Server (wcmd) as g:
    f = Memoize (f, DirStore (salt="water Dec 9"))
    g = Memoize (g, DirStore (salt=wcmd + "Dec 9"))

    # One could compose  (f + g, trafo), but we  use f(s) and g(s)
    # below as functions of internal coordinates:
    f = compose (f, trafo)
    g = compose (g, trafo)

    e = f + g

    print "s(start)=", s
    print "x(start)=\n", trafo (s)

    s, info = minimize (f, s, algo=1)
    print s, e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    s0 = s

    s, info = minimize (e, s, algo=1)
    print s, e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    s1 = s

    # SPC geometry
    watoms = read ("spc.xyz")
    x = watoms.get_positions ()
    s2 = trafo.pinv (x)

    ss = [s0, s1, s2]
    for s in ss:
        print "f=", f (s) / kcal, "g=", g (s) / kcal, e (s) / kcal, "(kcal)"

