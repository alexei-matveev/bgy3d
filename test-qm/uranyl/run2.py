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
from ase.io import read, write
from pts.memoize import Memoize, DirStore
from pts.units import kcal
from pts.func import compose, NumDiff
from pts.zmat import Fixed, Rigid, ManyBody
from pts.qfunc import QFunc
from pts.fopt import minimize
from rism import Server
from numpy import max, abs

cmd = \
"""
/home/alexei/darcs/bgy3d/guile/runbgy.scm
--solute "UO2_5H2O, SPC"
"""
#--norm-tol 1e-14 --dielectric 78.4
#--rho 0.0333295 --beta 1.6889 --L 160 --N 4096
#--solute "uranyl, SPC"
#--solvent "water, SPC/E"

calc = ParaGauss (cmdline="mpirun --np 4 ~/darcs/ttfs-mac/runqm",
                  input="uranyl+water.scm")

atoms = read ("uo22+,5h2o.xyz")

x = atoms.get_positions ()

xs = [x[3 * i: 3 * i + 3] for i in range (len (x) / 3)]
fs = [Rigid (y) for y in xs]
trafo = ManyBody (Fixed (xs[0]), *fs[1:])
s = trafo.pinv (x) # s == zeros (6 * len (fs) - 6)

assert max (abs (s)) == 0.0
assert max (abs (trafo (s) - x)) < 1.0e-10

with QFunc (atoms, calc) as f:
    with Server (cmd) as g:
        f = Memoize (f, DirStore (salt="uo22+, 5h2o, qm"))
        f = compose (f, trafo)
        g1 = Memoize (g, DirStore (salt="uo22+, rism"))
        g = lambda x: g1(x) * kcal
        g = compose (g, trafo)
        g = NumDiff (g)

        e = f # + g

        print "s(start)=", s
        print "x(start)=\n", trafo (s)
        print "e(start)=", e (s)  / kcal

        s, info = minimize (e, s, algo=1, maxit=50, ftol=1.0e-2, xtol=1.0e-2)

        print "converged =", info["converged"], "in", info["iterations"], "iterations"
        print "s(min)=", s
        print "x(min)=\n", trafo (s)
        atoms.set_positions (trafo (s))
        write ("a50.xyz", atoms)

        # units = [(kcal, "kcal"), (eV, "eV"), (Hartree, "Hartree")]
        # for u, uu in units:
        #     print "e = ", f(s)/u, "+", g(s)/u, "=", e(s)/u, "(%s)" % uu
        # for u, uu in units:
        #     print "g = ", f.fprime(s)/u, "+", g.fprime(s)/u, "=", e.fprime(s)/u, "(%s/Unit)" % uu
