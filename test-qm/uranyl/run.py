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

cmd = \
"""
/home/alexei/darcs/bgy3d/guile/runbgy.scm
--solvent "water, SPC/E" --rho 0.0333295 --beta 1.6889 --L 160 --N 4096
--solute "uranyl, SPC" --norm-tol 1e-14 --dielectric 78.4
"""

calc = ParaGauss(cmdline="mpirun -np 4 ~/darcs/ttfs-mac/runqm",
                       input="uranyl.scm")

atoms = read ("uo22+,gp.xyz")

x = atoms.get_positions ()

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

with QFunc (atoms, calc) as f:
    with Server (cmd) as g:
        f = Memoize (f, DirStore (salt="uo22+, gp, qm"))
        g = Memoize (g, DirStore (salt=cmd))

        # One could compose  (f + g, trafo), but we  use f(s) and g(s)
        # below as functions of internal coordinates:
        f = compose (f, trafo)
        g = compose (g, trafo)

        e = f + g

        print "s(start)=", s
        print "x(start)=\n", trafo (s)

        s, info = minimize (e, s, algo=1)

        print "converged =", info["converged"], "in", info["iterations"], "iterations"
        print "s(min)=", s
        print "x(min)=\n", trafo (s)

        units = [(kcal, "kcal"), (eV, "eV"), (Hartree, "Hartree")]
        for u, uu in units:
            print "e = ", f(s)/u, "+", g(s)/u, "=", e(s)/u, "(%s)" % uu
        for u, uu in units:
            print "g = ", f.fprime(s)/u, "+", g.fprime(s)/u, "=", e.fprime(s)/u, "(%s/Unit)" % uu
