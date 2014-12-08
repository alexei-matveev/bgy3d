#!/usr/bin/python -u

from __future__ import with_statement, print_function
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
from ase.io import read, write
from pts.memoize import Memoize, DirStore
from pts.units import kcal
from pts.func import compose
from pts.zmat import Fixed, Rigid, ManyBody #, Move, relate
from pts.cfunc import Cartesian
from pts.rc import Distance, Difference, Array
from pts.fopt import minimize, cminimize
from pts.path import Path
from rism import Server
from numpy import max, abs, zeros, array, asarray, linspace, savetxt, loadtxt, pi
print ("kcal=", kcal)

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

# With  this command  line the  RISM  server will  return the  "excess
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
--rmax=40
--nrad=1536
--closure=KH
--hnc
""" % (solvent_name, solute_name (NW))

atoms = read ("%dw,mm.xyz" % NW)

x = atoms.get_positions ()

# FIXME:  it  is assumed  here  that  all  species have  three  atoms,
# ordered:
xs = [x[3 * i: 3 * i + 3] for i in range (len (x) / 3)]
assert (len (xs) == 1 + NW)

def make_uranyl (x, flexible=False):
    """
    Return the "best" linear approximation to the uranyl geometry as a
    Fixed() func.
    """
    from pts.zmat import ZMat

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

    # Uranyl  may or  may not  be  fixed, but  to adjust  orientation,
    # rotate a rigid object:
    Y = Rigid (zmt (s))
    y = Y (Y.pinv (x))

    if not flexible:
        return Fixed (y)
    else:
        # return Move (relate (x, zmt (s)), zmt)
        return ManyBody (Fixed (x[0:1]), Cartesian (x[1:3]), dof=[0, 6])

flexible = True

if True:
    uranyl = make_uranyl (xs[0], flexible)
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

# FIXME: waters are free-floating, so formally for a flexible uranyl 3
# dofs  should  suffice. But  making  waters  dance  around uranyl  is
# counter-intuitive. Instead we set the OU oxygens free:
if flexible:
    dof = [6] + [6] * len (waters)
else:
    dof = [0] + [6] * len (waters)

trafo = ManyBody (uranyl, *waters, dof=dof)

# 6 dof per rigid water
if flexible:
    s = zeros (6 + 6 * len (waters))
    # s[0:3] = (1.79, 1.79, pi)
else:
    s = zeros (6 * len (waters))

# Destructively updates "atoms"
def write_xyz (path, x):
    atoms.set_positions (x)
    write (path, atoms)

#
# Bending force constant:
#
# Ref. [1]: 198 (1255) kJ/rad^2 = 2.05 (13.01) eV/rad^2
# Ref. [2]: 13.009 eV/rad^2
#
# Angle is 180, naturally. For  them kcal is 4.184 J. Stretching force
# constant:
#
# Ref. [1]: r0 = 0.176 nm, k = 622300 kJ/nm^2 = 64.50 eV/A^2
# Ref. [2]: r0 = 1.8 A, k = 43.364 eV/A^2
#
# [1] Tiwari et al., PCCP 16 8060 (2014)
#
# [2] Derived from GW force field in Kerisit and Liu,
#     J. Phys. Chem. A, 2013, 117 (30), pp 6421-6432
#
# Build  a  Func  of  cartesian coordinates,  stretching  and  bending
# parameters as tuples.  Parameters not hashed, beware of memoization!
#
if flexible:
    from pts.pes.ab2 import AB2
    # Base units here: A, radians, eV:
    if True:
        # PM13 soft uranyl:
        uo2 = AB2 ((1.76, 64.50), (pi, 2.05)) # Ref. [1]
    else:
        # KL13 aka GW hard uranyl:
        uo2 = AB2 ((1.80, 43.364), (pi, 13.009)) # Ref. [2]
else:
    uo2 = None


def test_uranyl ():
    e = compose (uo2, uranyl)
    # s = array ([1.79, 1.79, pi])
    s = zeros (6) + 0.1
    # print (uranyl (s))
    # exit (0)
    sm, info = minimize (e, s)
    print ("e(0)=", e (s), "e(1)=", e(sm), "iterations=", info["iterations"])

if uo2 is not None:
    test_uranyl()


def optimization (s):
    with Server (cmd) as g, Server (alt) as h:
        # Change  the  salts  to   discard  memoized  results  (or  delete
        # ./cache.d):
        g = Memoize (g, DirStore (salt=cmd + "TESTING1 Sep14"))
        h = Memoize (h, DirStore (salt=alt + "TESTING1 Sep14"))

        # A  Func  of cartesian  coordinates,  stretching and  bending
        # parameters  as  tuples.  Parameters  not  hashed, beware  of
        # memoization of this term!:
        if uo2 is not None:
            g = g + uo2

        g = compose (g, trafo)  # MM self-energy
        h = compose (h, trafo)  # RISM solvation energy

        def opt (e, s, name, **kwargs):
            print ("XXX: " + name + "...")

            print (name + ": e(0)=", e (s) / kcal, "kcal, |g|=", max (abs (e.fprime (s))) / kcal, "kcal/Unit")

            s, info = minimize (e, s, **kwargs)

            print (name + ": converged =", info["converged"], "in", info["iterations"], "iterations")
            write_xyz (name + ".xyz", trafo (s))

            # print info
            if "trajectory" in info:
                traj = info["trajectory"]
                print (map(e, traj))
                for i, si in enumerate (traj):
                    # write_xyz ("traj-%s-%03d.xyz" % (name, i), trafo (si))
                    pass

            return s, info


        # MM self-energy:
        with g as e:
            s0, info = opt (e, s, "MM", algo=1, maxstep=0.1, maxit=200, ftol=5.0e-3, xtol=5.0e-3)
            print ("XXX: MM", e (s), "eV", e (s) / kcal, "kcal", info["converged"])

        exit (0)
        # MM self-energy with RISM solvation:
        with g + h as e:
            s1, info = opt (e, (s if NW == 4 else s0), "MM+RISM", algo=1, maxstep=0.1, maxit=200, ftol=5.0e-3, xtol=5.0e-3)
            print ("XXX: MM+RISM", e (s), "eV", e (s) / kcal, "kcal", info["converged"])

        # This prints  a table of  various functionals applied  to several
        # geometries:
        ss = [s0, s1]
        with g + h as e:
            for s in ss:
                print ("MM=", g (s) / kcal, "RISM=", h (s) / kcal, "MM+RISM=", e (s) / kcal, "(kcal)")


# optimization (s)
# exit (0)

def initial_path (f, s, c):
    """
    This  is  called  with  pure  MM  PES f(s)  to  get  some  initial
    path. Somewhat ad-hoc.
    """
    # The input  geometry s may  be reasonable, but  it may not  be an
    # exact minimum. Reoptimize:
    s, info = minimize (f, s, maxit=100, ftol=5.0e-4, xtol=5.0e-4, algo=1)
    print ("converged=", info["converged"], "in", info["iterations"])

    # The value of reaction coordinate in the initial geometry:
    c0 = c(s)
    print ("XXX: rc(0)=", c0)

    it = 0
    ss = []
    qs = []
    algo = 1
    while c(s) > -c0:
        it += 1
        # write_xyz ("in-%03d.xyz" % it, trafo (s))
        print ("XXX: rc(0)=", c(s))
        s, info = cminimize (f, s, Array (c), maxit=200, ftol=5.0e-3, xtol=5.0e-3, algo=algo)
        # Collect optimized geometries and the corresponding
        # values of reaction coordinate:
        ss.append (s)
        qs.append (c(s))
        print ("converged=", info["converged"], "in", info["iterations"])
        print ("XXX: rc(1)=", c(s))
        # write_xyz ("out-%03d.xyz" % it, trafo (s))
        # Here  c.fprime(s) is how  much that  reaction coordinate
        # will change if you  modify the coordinates. Hacky way to
        # chage a geometry so that the RC is modified too:
        s = s - 0.1 * c.fprime (s)

    ss = asarray (ss)
    qs = asarray (qs)
    print ("qs=", qs)
    p = Path (ss, qs)
    print ("path=", p)
    qs = linspace (c0, -c0, 21)
    print ("qs=", qs)
    ss = map (p, qs)
    res = map (lambda s: cminimize (f, s, Array (c), maxit=200, ftol=1.0e-3, xtol=1.0e-3, algo=algo), ss)
    infos = [inf for _, inf in res]
    ss = [s for s, _ in res]
    for info in infos:
        print ("converged=", info["converged"], "in", info["iterations"])
    for i, s in enumerate (ss):
        write_xyz ("out-%03d.xyz" % i, trafo (s))

    return qs, ss


def exchange (s):
    # Functions of cartesian coordinates:
    ra = Distance ([0, 3])
    rb = Distance ([0, 18])
    rc = Difference (ra, rb)
    with Server (cmd) as f0, Server (alt) as h0:
        f1 = Memoize (f0, DirStore (salt=cmd + "TESTING"))
        h1 = Memoize (h0, DirStore (salt=alt + "TESTING"))

        # A  Func  of cartesian  coordinates,  stretching and  bending
        # parameters  as  tuples.  Parameters  not  hashed, beware  of
        # memoization of this term!:
        if uo2 is not None:
            f1 = f1 + uo2

        f = compose (f1, trafo)
        h = compose (h1, trafo)
        c = compose (rc, trafo)

        refine = True
        if refine:
            ss = loadtxt ("ss,initial.txt")
            qs = array (map (c, ss))
        else:
            qs, ss = initial_path(f, s, c)

        p = Path (ss, qs)

        if True:
            # Energy profile, smooth:
            print ("# q, ra(q), rb(q), E(q)")
            for q in linspace (qs[0], qs[-1], 100):
                print (q, ra (trafo (p(q))), rb (trafo (p (q))), f(p(q)))

        if True:
            # Energy profile, coarse, with dG(q):
            with open ("./profile,intial.txt", "w") as file:
                print ("# q, ra(q), rb(q), E(q), dG(q)", file=file)
                for q in qs:
                    print (q, ra (trafo (p(q))), rb (trafo (p (q))), f(p(q)), h(p(q)), file=file)

        with f + h as e:
            # Terminals  were constrained, obtain  fully unconstrained
            # ones:
            if refine:
                sab = loadtxt ("sab,initial.txt")
            else:
                sab = (ss[0], ss[-1])

            def opt (s):
                sm, info = minimize (e, s, maxit=50, ftol=1.0e-2, xtol=1.0e-2, algo=1)
                print ("converged=", info["converged"], "in", info["iterations"])
                return sm
            sab = map (opt, sab)
            qab = map (c, sab)
            qa, qb = qab
            sa, sb = sab
            savetxt ("sab.txt", sab)
            write_xyz ("KH-aaa.xyz", trafo (sa))
            write_xyz ("KH-bbb.xyz", trafo (sb))

            def copt (s):
                sm, info = cminimize (e, s, Array (c), maxit=50, ftol=1.0e-2, xtol=1.0e-2, algo=0)
                print ("converged=", info["converged"], "in", info["iterations"])
                return sm

            ss = array (map (copt, ss))
            savetxt ("ss.txt", ss)

            # Energy profile:
            with open ("./profile,rs.txt", "w") as file:
                print ("# Optimized profile:", file=file)
                print (qa, f(sa), h(sa), e(sa), file=file)
                for q, s in zip (qs, ss):
                    print (q, f(s), h(s), e(s), file=file)
                print (qb, f(sb), h(sb), e(sb), file=file)

            for i, s in enumerate (ss):
                write_xyz ("KH-%03d.xyz" % i, trafo (s))

            p = Path (ss, qs)
            with open ("./profile,rs,interp.txt", "w") as file:
                print ("# Interpolated profile:", file=file)
                for q in linspace (qs[0], qs[-1], 3 * len (qs)):
                    print (q, f(p(q)), h(p(q)), e(p(q)), file=file)


exchange (s)


