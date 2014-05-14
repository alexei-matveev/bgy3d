#!/usr/bin/python

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
# Note that this driver script communicates with MPI processes using a
# named FIFO.  For that to work  the driver script must be executed on
# the same node as the rank-0 MPI process. This is NOT guaranteed when
# one  uses "salloc  -n8 run.py".   It  is an  aknowledged feature  of
# "salloc" that the  driver script is executed locally  and only those
# commands  prefixed  by  srun/mpirun  use the  allocated  nodes.  Use
# "sbatch -n8 run.py" instead.
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
mpirun /users/alexei/darcs/bgy3d-wheezy/guile/runbgy.scm
--solute "UO2_5H2O, SPC"
"""

# With  thin command  line the  RISM  server will  return the  "excess
# chemical  potential" of  the  solute in  solvent.  This is,  roughly
# speaking, the "averaged" solute-solvent interaction:
alt = \
"""
mpirun /users/alexei/darcs/bgy3d-wheezy/guile/runbgy.scm
--solvent "water, SPC/E"
--solute "UO2_5H2O, SPC"
--norm-tol 1e-14 --dielectric 78.4
--rho 0.0333295 --beta 1.6889 --L 20 --N 512
"""
alt1 = \
"""
mpirun /users/alexei/darcs/bgy3d-wheezy/guile/runbgy.scm
--solvent "water, SPC/E"
--solute "UO2_5H2O, SPC"
--norm-tol 1e-14 --dielectric 78.4
--rho 0.0333295 --beta 1.6889 --L 160 --N 4096
"""

calc = ParaGauss (cmdline="mpirun ~/darcs/ttfs-mac/runqm",
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

def read_xyz (path):
    atoms = read (path)
    return atoms.get_positions ()

def make_waters (xs):
    # All waters will be rigid objects with the same geometry:
    w = Rigid (read_xyz ("spc.xyz"))
    # print xs
    # print read_xyz ("spc.xyz")

    # But at different locations and orientation:
    ss = [w.pinv (x) for x in xs]

    return [Rigid (w (s)) for s in ss]


waters = make_waters (xs[1:])
# waters = [Rigid (y) for y in xs[1:]]

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

# Destructively updates "atoms"
def write_xyz (path, x):
    atoms.set_positions (x)
    write (path, atoms)

with QFunc (atoms, calc) as f, Server (cmd) as g, Server (alt) as h, Server (alt1) as h1:
    # Change  the  salts  to   discard  memoized  results  (or  delete
    # ./cache.d):
    f = Memoize (f, DirStore (salt="uo22+, 5h2o, bp"))
    g = Memoize (g, DirStore (salt=cmd + "Dec3(b)"))
    h = Memoize (h, DirStore (salt=alt + "Dec3(b)"))
    h1 = Memoize (h1, DirStore (salt=alt1 + "Dec 17a"))

    f = compose (f, trafo)  # QM self-energy
    g = compose (g, trafo)  # MM self-energy
    h = compose (h, trafo)  # RISM solvation energy
    h1 = compose (h1, trafo)  # RISM solvation energy, hi precision

    def opt (e, s, name, **kwargs):
        print "XXX: " + name + "..."

        print name + ": e(0)=", e (s) / kcal, "kcal, |g|=", max (abs (e.fprime (s))) / kcal, "kcal/Unit"

        s, info = minimize (e, s, **kwargs)

        print name + ": converged =", info["converged"], "in", info["iterations"], "iterations"
        write_xyz (name + ".xyz", trafo (s))

        # print info
        traj = info["trajectory"]
        print map(e, traj)
        for i, si in enumerate (traj):
            # write_xyz ("traj-%s-%03d.xyz" % (name, i), trafo (si))
            pass

        return s, info


    # MM self-energy:
    with g as e:
        s, info = opt (e, s, "MM", algo=1, maxstep=0.1, maxit=100, ftol=5.0e-3, xtol=5.0e-3)
        print "XXX: MM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    s0 = s

    # MM self-energy with RISM solvation:
    with g + h as e:
        s, info = opt (e, s, "MM+RISM", algo=1, maxstep=0.1, maxit=100, ftol=5.0e-3, xtol=5.0e-3)
        print "XXX: MM+RISM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    s1 = s

    # QM self-energy (start with MM geom):
    with f as e:
        s, info = opt (e, s0, "QM", algo=1, maxstep=0.1, maxit=100, ftol=1.0e-2, xtol=1.0e-2)
        print "XXX: QM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    s2 = s

    # QM self-energy with RISM solvation (diverges):
    with f + h as e:
        s, info = opt (e, s1, "QM+RISM", algo=1, maxstep=0.1, maxit=98, ftol=5.0e-3, xtol=5.0e-3)
        print "XXX: QM+RISM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    s3 = s

    ss = [s0, s1, s2, s3]
    with f + h as e, g + h as e1:
        for s in ss:
            print "QM=", f (s) / kcal, "MM=", g (s) / kcal, "RISM=", h (s) / kcal, \
                "QM+RISM=", e (s) / kcal, "MM+RISM=", e1 (s) / kcal, "(kcal)"


    clean ()

    # MM self-energy with RISM solvation:
    with g + h1 as e:
        s, info = opt (e, s1, "uranyl+water,mm+rism", algo=1, maxstep=0.1, maxit=100, ftol=1.0e-2, xtol=1.0e-2)
        print "XXX: MM+RISM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    S1 = s

    # QM self-energy with RISM solvation (diverges):
    # with f + h1 as e:
    #     s, info = opt (e, s3, "uranyl+water,qm+rism", algo=1, maxstep=0.1, maxit=100, ftol=1.0e-2, xtol=1.0e-2)
    #     print "XXX: QM+RISM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]
    # S3 = s

    ss = [s0, S1, s2, s3]
    with h1 as h:
        with f + h as e, g + h as e1:
            for s in ss:
                print "QM=", f (s) / kcal, "MM=", g (s) / kcal, "RISM=", h (s) / kcal, \
                    "QM+RISM=", e (s) / kcal, "MM+RISM=", e1 (s) / kcal, "(kcal)"
