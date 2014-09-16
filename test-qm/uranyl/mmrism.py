#!/usr/bin/python

from __future__ import with_statement
#
# Tell Python interpreter where to find custom modules:
#
#   export PYTHONPATH=~/PYTHON:~/darcs/bgy3d-wheezy/python
#
# For the runbgy.scm  script to run one might need  to set the library
# search path:
#
#   export LD_LIBRARY_PATH=~/darcs/bgy3d-wheezy
#
# Note that this driver script communicates with MPI processes using a
# named FIFO.  For that to work  the driver script must be executed on
# the same node as the rank-0 MPI process. This is NOT guaranteed when
# one uses "salloc -n8 ./mmrism.py".   It is an aknowledged feature of
# "salloc" that the  driver script is executed locally  and only those
# commands prefixed by srun/mpiexec  use the allocated nodes.  Instead
# use
#
#   sbatch --exclude=quad4 -n 8 ./mmrism.py
#
# This script is not self-contained, it refers to other files that are
# used as input:
#
#   uo22+,%dh2o.xyz
#
#     Initial geometry. It does  not necessarily use the same geometry
#     for all water molecules. The uranyl may be slightly bent.
#
#   spc.xyz
#
#     Geometry of SPC/E water.
#
import os
from ase.io import read, write
from pts.memoize import Memoize, DirStore
from pts.units import kcal, eV, Hartree
from pts.func import compose
from pts.zmat import Fixed, Rigid, ManyBody
from pts.fopt import minimize
from rism import Server
from numpy import max, abs, zeros, array

# File  names  and  command  line   flags  depend  on  the  number  of
# waters. FIXME: name solutes uniformely:
NW = 6
def solute_name (nw):
    if nw == 4:
        return "UO2_%dH2O, D4H, KL2-PRSPC" % nw
    elif nw == 6:
        return "UO2_%dH2O, D3D, KL2-PRSPC" % nw
    else:
        return "UO2_%dH2O, KL2-PRSPC" % nw

solvent_name = "water, PR-SPC/E"

# With this command  line the RISM server will  return the "gas phase"
# or self-energy of the solute:
cmd = \
"""
mpiexec /users/alexei/darcs/bgy3d-wheezy/guile/runbgy.scm
--solute "%s"
""" % solute_name (NW)

# With  thin command  line the  RISM  server will  return the  "excess
# chemical  potential" of  the  solute in  solvent.  This is,  roughly
# speaking, the "averaged" solute-solvent interaction:
alt = \
"""
mpiexec /users/alexei/darcs/bgy3d-wheezy/guile/runbgy.scm
--solvent "%s"
--solute "%s"
--norm-tol=1e-14
--dielectric=78.4
--rho=0.0333295
--beta=1.6889
--L=10
--N=96
--closure=KH
--hnc
""" % (solvent_name, solute_name (NW))

atoms = read ("uo22+,%dh2o.xyz" % NW)

x = atoms.get_positions ()

# FIXME:  it  is assumed  here  that  all  species have  three  atoms,
# ordered:
xs = [x[3 * i: 3 * i + 3] for i in range (len (x) / 3)]
assert (len (xs) == 1 + NW)

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

    # Set the bond  lengths equal and 180 degrees  angle. Note that if
    # you plug these internal variables into z-matrix you will not get
    # zmt(s) ~  x even when s was  derived from x. That  is because of
    # the default orientation z-matrix chooses:
    s_bond = 1.79 # A or sum (s[:2]) / 2
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

# Destructively updates "atoms"
def write_xyz (path, x):
    atoms.set_positions (x)
    write (path, atoms)

with Server (cmd) as g, Server (alt) as h:
    # Change  the  salts  to   discard  memoized  results  (or  delete
    # ./cache.d):
    g = Memoize (g, DirStore (salt=cmd + "Sep14"))
    h = Memoize (h, DirStore (salt=alt + "Sep14"))

    g = compose (g, trafo)  # MM self-energy
    h = compose (h, trafo)  # RISM solvation energy

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
        s0, info = opt (e, s, "MM", algo=1, maxstep=0.1, maxit=200, ftol=5.0e-3, xtol=5.0e-3)
        print "XXX: MM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]

    # MM self-energy with RISM solvation:
    with g + h as e:
        s1, info = opt (e, (s if NW == 4 else s0), "MM+RISM", algo=1, maxstep=0.1, maxit=200, ftol=5.0e-3, xtol=5.0e-3)
        print "XXX: MM+RISM", e (s), "eV", e (s) / kcal, "kcal", info["converged"]

    # This prints  a table of  various functionals applied  to several
    # geometries:
    ss = [s0, s1]
    with g + h as e:
        for s in ss:
            print "MM=", g (s) / kcal, "RISM=", h (s) / kcal, "MM+RISM=", e (s) / kcal, "(kcal)"
    exit (0)
