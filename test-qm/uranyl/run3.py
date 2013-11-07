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
from pts.units import kcal, Hartree, eV
from pts.func import compose, NumDiff
from pts.cfunc import Affine
from pts.zmat import ZMat, Rigid, Fixed, ManyBody, Move, relate
from pts.qfunc import QFunc
from pts.fopt import minimize
from rism import Server
from numpy import max, abs, zeros, vstack, pi, empty

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

command = "salloc -n16 mpirun ~/darcs/ttfs-mac/runqm"

calc = ParaGauss (cmdline=command, input="water.scm")

atoms = read ("h2o.xyz")

x = atoms.get_positions ()

#
# Z-matrix for water:
#
zmt = ZMat ([(None, None, None),
             (1, None, None),
             (1, 2, None)], base=1)

#
# Initial values of internal coordinates:
#
s = zmt.pinv (x)
assert max (abs (s - zmt.pinv (zmt (s)))) < 1.0e-10

clean ()

with QFunc (atoms, calc) as f:
    f = Memoize (f, DirStore (salt="h2o, qm"))
    e = compose (f, zmt)

    print s, e (s)
    s, info = minimize (e, s, algo=1, ftol=1.0e-2, xtol=1.0e-2)
    print s, e (s), info["converged"]

#
# Internal coordinates:
#
# In:  0.97843364  0.97826543  1.87152024
# Out: 0.98477091  0.98477695  1.88111918
#
# Rigid water with optimized internal coordinates:
#
X = Rigid (zmt (s))

# Five optimized  rigid waters, actual location and  orientation to be
# specified ...
waters = [X] * 5

atoms = read ("uo22+,5h2o.xyz")
xa = atoms.get_positions ()

# Bond length and the bond angle for uranyl, use the same Z-matrix for
# uranyl:
s = zmt.pinv (xa[:3])

# Average bond length and 180 degrees angle:
s_bond = sum (s[:2]) / 2
s = [s_bond, s_bond, pi]

# Uranyl  will be  fixed, but  to adjust  orientation, rotate  a rigid
# object:
Y = Rigid (zmt (s))
uranyl = Fixed (Y (Y.pinv (xa[:3])))

trafo = ManyBody (uranyl, *waters)

# Here  initial  positions  and  orientation  of 5  rigid  waters  are
# computed to  approximate the geometry  from the xyz-file.  Uranyl is
# fixed:
s = trafo.pinv (xa)

calc = ParaGauss (cmdline=command, input="uranyl+water.scm")

clean ()

with QFunc (atoms, calc) as f:
    f = Memoize (f, DirStore (salt="uranyl+water, qm"))
    e = compose (f, trafo)

    print "s(start)=", s
    print "x(start)=\n", trafo (s)
    print "e(start)=", e (s)  / kcal

    s, info = minimize (e, s, algo=1, ftol=1.0e-2, xtol=1.0e-2)

    print "converged =", info["converged"], "in", info["iterations"], "iterations"
    print "s(min)=", s
    print "x(min)=\n", trafo (s)
    atoms.set_positions (trafo (s))
    write ("a50.xyz", atoms)

    units = [(kcal, "kcal"), (eV, "eV"), (Hartree, "Hartree")]
    for u, uu in units:
        print "e = ", e(s) / u, "%s," % uu, \
            "|g| = ", max(abs(e.fprime(s))) / u, "%s/Unit" % uu


xa = trafo (s)
print "x(min=)\n", xa

# Bond length and the bond angle for uranyl, use the same Z-matrix for
# uranyl:
s = zmt.pinv (xa[:3])
print "s(min)=", s

# Average bond length:
s_bond = sum (s[:2]) / 2
s = [s_bond]

# Geometry  of uranyl  as a  function  of single  argument. FIXME:  no
# special case for scalar argument!
def uranyl (v):
    z, = v
    return [[0, 0, 0],
            [0, 0, z],
            [0, 0, -z]]

# A linear Func() of a 1-array:
uranyl = Affine (uranyl, s)

# Relative position of second geometry wrt the first as a 6-array:
q = relate (uranyl (s), xa[:3])

# Move/rotate output of uranyl() Func():
uranyl = Move (q, uranyl)

# Flexible uranyl and five rigid waters at their (more or less) proper
# positions.   Each water is  characterised by  6 degrees  of freedom,
# uranyl by just one:
trafo = ManyBody (uranyl, *waters, dof=[1, 6, 6, 6, 6, 6])

# Initial values for all of degrees of freedom:
s = trafo.pinv (xa)
print "diff\n", max (abs (trafo (s) - xa))

clean ()

with QFunc (atoms, calc) as f:
    f = Memoize (f, DirStore (salt="uranyl+water, qm"))
    e = compose (f, trafo)

    print "s(start)=", s
    print "x(start)=\n", trafo (s)
    print "e(start)=", e (s)  / kcal

    s, info = minimize (e, s, algo=1, ftol=1.0e-2, xtol=1.0e-2)

    print "converged =", info["converged"], "in", info["iterations"], "iterations"
    print "s(min)=", s
    print "x(min)=\n", trafo (s)
    atoms.set_positions (trafo (s))
    write ("a50.xyz", atoms)

    units = [(kcal, "kcal"), (eV, "eV"), (Hartree, "Hartree")]
    for u, uu in units:
        print "e = ", e(s) / u, "%s," % uu, \
            "|g| = ", max(abs(e.fprime(s))) / u, "%s/Unit" % uu
