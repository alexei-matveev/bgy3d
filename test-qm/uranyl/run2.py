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
from ase.io import read, write
from pts.memoize import Memoize, DirStore
from pts.units import kcal, eV, Hartree
from pts.func import compose, NumDiff
from pts.zmat import Fixed, Rigid, ManyBody
from pts.qfunc import QFunc
from pts.fopt import minimize
from rism import Server
from numpy import max, abs, zeros

# With this command  line the RISM server will  return the "gas phase"
# or self-energy of the solute:
cmd = \
"""
/home/alexei/darcs/bgy3d/guile/runbgy.scm
--solute "UO2_5H2O, SPC"
"""

# With  thin command  line the  RISM  server will  return the  "excess
# chemical  potential" of  the  solute in  solvent.  This is,  roughly
# speaking, the "averaged" solute-solvent interaction:
alt = \
"""
/home/alexei/darcs/bgy3d/guile/runbgy.scm
--solvent "water, SPC/E"
--solute "UO2_5H2O, SPC"
--norm-tol 1e-14 --dielectric 78.4
--rho 0.0333295 --beta 1.6889 --L 20 --N 512
"""

calc = ParaGauss (cmdline="mpirun --np 4 ~/darcs/ttfs-mac/runqm",
                  input="uranyl+water.scm")

atoms = read ("uo22+,5h2o.xyz")
# atoms = read ("a50.xyz")

x = atoms.get_positions ()

xs = [x[3 * i: 3 * i + 3] for i in range (len (x) / 3)]

def make_uranyl (x):
    """
    Return the "best" linear approximation to the uranyl geometry as a
    Fixed() func.
    """
    from pts.zmat import ZMat
    from numpy import pi

    # Make uranyl linear:
    zmt = ZMat ([(None, None, None),
                 (1, None, None),
                 (1, 2, None)], base=1)

    # Bond length and the bond angle for uranyl, use the same Z-matrix
    # for uranyl:
    s = zmt.pinv (x)

    # Average bond length and 180 degrees angle. Note that if you plug
    # these internal variables into z-matrix you will not get zmt(s) ~
    # x  even though  s was  derived from  x. That  is because  of the
    # default orientation z-matrix chooses:
    s_bond = sum (s[:2]) / 2
    s = [s_bond, s_bond, pi]

    # Uranyl  will be  fixed, but  to adjust  orientation, rotate  a rigid
    # object:
    Y = Rigid (zmt (s))
    return Fixed (Y (Y.pinv (x)))

if True:
    uranyl = make_uranyl (xs[0]) # linear
else:
    uranyl = Fixed (xs[0])      # may be bent

waters = [Rigid (y) for y in xs[1:]]
trafo = ManyBody (uranyl, *waters)

# 6 dof per rigid water:
s = zeros (6 * len (waters))

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

clean ()


with QFunc (atoms, calc) as f, Server (cmd) as g, Server (alt) as h:
    # Change  the  salts  to   discard  memoized  results  (or  delete
    # ./cache.d):
    f = Memoize (f, DirStore (salt="uo22+, 5h2o, qm"))
    g = Memoize (g, DirStore (salt=cmd + "Dec3(b)"))
    h = Memoize (h, DirStore (salt=alt + "Dec3(b)"))

    f = compose (f, trafo)  # QM self-energy
    g = compose (g, trafo)  # MM self-energy
    h = compose (h, trafo)  # RISM solvation energy

    e = g + h # + f

    print "s(start)=", s
    print "x(start)=\n", trafo (s)
    print "e(start)=", e (s) / kcal, "g(start)=", g (s) / kcal, "h(start)=", h (s) / kcal

    s, info = minimize (e, s, algo=1, maxstep=0.1, maxit=100, ftol=1.0e-2, xtol=1.0e-2)

    print "converged =", info["converged"], "in", info["iterations"], "iterations"
    print "s(min)=", s
    print "x(min)=\n", trafo (s)
    atoms.set_positions (trafo (s))
    write ("a50.xyz", atoms)

    units = [(kcal, "kcal"), (eV, "eV"), (Hartree, "Hartree")]
    for u, uu in units:
        print "e = ", g(s)/u, "+", h(s)/u, "=", e(s)/u, "(%s)" % uu

    for u, uu in units:
        print "e = ", e(s) / u, "%s," % uu, \
            "|g| = ", max(abs(e.fprime(s))) / u, "%s/Unit" % uu


    e = f + h
    for u, uu in units:
        print "e = ", e(s) / u, "%s," % uu, \
            "|g| = ", max(abs(e.fprime(s))) / u, "%s/Unit" % uu
