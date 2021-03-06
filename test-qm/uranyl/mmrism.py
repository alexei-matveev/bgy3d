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
from ase.calculators.paragauss import ParaGauss
from pts.memoize import Memoize, DirStore
from pts.units import kcal
from pts.func import compose, Inverse
from pts.qfunc import QFunc
from pts.zmat import Fixed, Rigid, ManyBody #, Move, relate
from pts.cfunc import Cartesian, Affine
from pts.rc import Distance, Difference, Array
from pts.fopt import minimize, cminimize
from pts.path import Path, MetricPath
from pts.metric import Metric
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
mpiexec /users/alexei/darcs/bgy3d/guile/runbgy.scm
--solute "%s"
""" % solute_name (NW)

# With  this command  line the  RISM  server will  return the  "excess
# chemical  potential" of  the  solute in  solvent.  This is,  roughly
# speaking, the "averaged" solute-solvent interaction:
alt = \
"""
mpiexec /users/alexei/darcs/bgy3d/guile/runbgy.scm
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

# Command  line for  ./sym/runqmmm  used as  PG  calculator. Input  is
# specified separately on call to constructor.
qmd0 = \
"""
mpiexec ./sym/runqmmm
--solvent "water, PR-SPC/E"
--solute "uranyl, 6w, pcm"
--norm-tol=1e-14
--dielectric=78.4
--rho=0.0333295
--beta=1.6889
--L=10
--N=64
--rmax=40
--nrad=1536
--closure=KH
--hnc
"""

qmd = "mpiexec /users/alexei/git/ttfs-work-gpl/runqm"

atoms = read ("%dw,mm.xyz" % NW)

x = atoms.get_positions ()

# FIXME:  it  is assumed  here  that  all  species have  three  atoms,
# ordered:
xs = [x[3 * i: 3 * i + 3] for i in range (len (x) / 3)]
assert (len (xs) == 1 + NW)

def make_uranyl (x, flexible=6):
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

    if flexible == 0:
        return Fixed (y)    # 0 dof
    elif flexible == 2:
        # liear uranyl with flexible bonds:
        def f (x):
            ra, rb = x
            return array ([[0., 0., 0.],
                           [0., 0., ra],
                           [0., 0., -rb]])
        return Affine (f, array ([0., 0.])) # 2 dof
    elif flexible == 6:
        # return Move (relate (x, zmt (s)), zmt)
        return ManyBody (Fixed (x[0:1]), Cartesian (x[1:3]), dof=[0, 6])
    else:
        assert False

flexible = 2

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
# counter-intuitive. Instead we set the OU oxygens free. For isotropic
# PES  this, however,  introduces  degeneracy wrt  orientation of  the
# complex as a whole. The path images, for example, may "diffuse" from
# each  other  in  orientation   so  that  interpolation  will  become
# problematic.
dof = [flexible] + [6] * len (waters)

trafo = ManyBody (uranyl, *waters, dof=dof)

# 6 dof per rigid water
s = zeros (flexible + 6 * len (waters))

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
if flexible != 0:
    from pts.pes.ab2 import AB2
    # Base units here: A, radians, eV:
    if True:
        # PM13 soft uranyl:
        uo2 = AB2 ((1.76, 64.50), (pi, 2.05)) # Ref. [1]
    if False:
        # KL13 aka GW hard uranyl:
        uo2 = AB2 ((1.80, 43.364), (pi, 13.009)) # Ref. [2]
    if False:
        print ("CMDLINE=", qmd)
        uo2 = QFunc (atoms, ParaGauss (cmdline=qmd, input="6w,c1.scm"))
else:
    uo2 = None


def test_uranyl ():
    e = compose (uo2, uranyl)
    e = Memoize (e, DirStore (salt="XXX YYY"))
    # s = array ([1.79, 1.79, pi])
    s = zeros (6) # + 0.1
    print (uranyl (s))
    # print (e(s))
    # exit (0)
    sm, info = minimize (e, s, ftol=1.0e-2, xtol=1.0e-2, algo=1)
    print (uranyl (sm))
    print ("e(0)=", e (s), "e(1)=", e(sm), "iterations=", info["iterations"])

if False:
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

def make_bias (r0, k):
    from pts.pes.bias import Bias
    return \
        Bias (r0, k , [0, 6]) + \
        Bias (r0, k , [0, 9]) + \
        Bias (r0, k , [0, 12]) + \
        Bias (r0, k , [0, 15]) + \
        Bias (r0, k , [0, 18])

#
# B(r) =  25 * (r - 2.50)^2  eV, for r >  2.50 A. For 2.51  A the bias
# amounts to  2.5 meV or 0.06  kcal. For 2.55  A it is 25  times more:
# 0.0625 eV or 1.44 kcal.
#
bias = make_bias (2.50, 50.0)

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
    # rc = Difference (ra, rb)
    rc = ra
    f0 = Server (cmd)
    h0 = Server (alt)
    F0 = QFunc (atoms, ParaGauss (cmdline=qmd, input="6w,c1.scm"))
    with f0 as f0, h0 as h0, F0 as F0:
        f1 = Memoize (f0, DirStore (salt=cmd + "TESTING"))
        h1 = Memoize (h0, DirStore (salt=alt + "TESTING"))
        F1 = Memoize (F0, DirStore (salt=qmd + "XXX4 50/291"))

        # A  Func  of cartesian  coordinates,  stretching and  bending
        # parameters  as  tuples.  Parameters  not  hashed, beware  of
        # memoization of this term!:
        if uo2 is not None:
            f1 = f1 + uo2

        f = compose (f1, trafo)
        h = compose (h1, trafo)
        c = compose (rc, trafo)
        F = compose (F1, trafo)
        B = compose (bias, trafo)

        refine = True
        if False:
            ss = loadtxt ("ss,initial.txt")
            qs = array (map (c, ss))
            print ("qs=", qs)
            print ("dq=", qs[1:] - qs[:-1])
            if False:
                P = Path (ss, qs)
                C = compose (c, P)
                print ("q=", qs)
                def solve (q):
                    x = q
                    while abs (C(x) - q) > 1.0e-10:
                        f, fprime = C.taylor (x)
                        df = f - q
                        x = x - df / fprime
                    return x
                xs = []
                for q in linspace (qs[0], qs[-1], 21):
                    x = solve (q)
                    print ("x=", x, "C(x)=", C(x), "q=", q)
                    xs.append (x)
                ss = array (map (P, xs))
                qs = array (map (c, ss))
                print ("qs=", qs)
                print ("dq=", qs[1:] - qs[:-1])
            p = Path (ss, qs)

        if False:
            qs, ss = initial_path(f + B, s, c)
            savetxt ("ss-initial.txt", ss)
            p = Path (ss, qs)

        if False:
            for i, q in enumerate (qs):
                write_xyz ("in-%03d.xyz" % i, trafo (p (q)))

        if False:
            # Energy profile, smooth:
            with open ("profile-initial.txt", "w") as file:
                print ("# q, ra(q), rb(q), E(q), E_qm(q), B(q)", file=file)
                for q in linspace (qs[0], qs[-1], 100):
                    print (q, ra (trafo (p(q))), rb (trafo (p (q))), f(p(q)), B (p(q)), file=file)

        if False:
            # Energy profile, coarse, with dG(q):
            with open ("./profile,intial.txt", "w") as file:
                print ("# q, ra(q), rb(q), F(q)", file=file)
                # print ("# q, ra(q), rb(q), E(q), dG(q), F(q)", file=file)
                for q in qs:
                    print (q, ra (trafo (p(q))), rb (trafo (p (q))), F(p(q)), file=file)
                    # print (q, ra (trafo (p(q))), rb (trafo (p (q))), f(p(q)), h(p(q)), F(p(q)), file=file)

        # with f + h as e:
        with F + h + B as e:
            print ("XXX")
            if False:
                p1 = Path (loadtxt ("sxy.txt", ndmin=2))
                c1 = compose (c, p1)
                q1 = linspace (0., 1., 21)
                for i, q in enumerate (q1):
                    write_xyz ("interp-%03d.xyz" % i, trafo (p1 (q)))
                print ("c=", [c1(q) for q in q1])
                Q1 = Inverse (c1)
                q2 = [Q1(c_) for c_ in [3.05]]
                print ("q=", q2, "c=", [c1(q) for q in q2])
                s2 = array ([p1(q) for q in q2])
                savetxt ("sxy1.txt", s2)
                exit (0)

            def opt (s):
                sm, info = minimize (e, s, maxit=50, ftol=1.0e-2, xtol=1.0e-2, algo=1)
                print ("converged=", info["converged"], "in", info["iterations"])
                return sm

            if True:
                # Unconstrained termnals:
                sab = loadtxt ("sab,initial.txt", ndmin=2)
                # sab = (ss[0], ss[-1])

                print ("before: q=", map (c, sab))
                sab = map (opt, sab)
                print ("after: q=", map (c, sab))
                savetxt ("sab.txt", sab)

                for i, s in enumerate (sab):
                    write_xyz ("min-%03d.xyz" % i, trafo (s))
                exit (0)

            def copt (s, maxit):
                sm, info = cminimize (e, s, Array (c), maxit=maxit, ftol=1.0e-3, xtol=1.0e-3, algo=0)
                print ("converged=", info["converged"], "in", info["iterations"])
                return sm

            maxit = [30] * len (ss)
            ss = array (map (copt, ss, maxit))
            savetxt ("ss.txt", ss)

            # Energy profile:
            with open ("./profile,rs.txt", "w") as file:
                print ("# Optimized profile:", file=file)
                # print (qa, f(sa), h(sa), e(sa), file=file)
                for q, s in zip (qs, ss):
                    print (q, f(s), h(s), e(s), file=file)
                # print (qb, f(sb), h(sb), e(sb), file=file)

            for i, s in enumerate (ss):
                write_xyz ("KH-%03d.xyz" % i, trafo (s))

            if False:
                p = Path (ss, qs)
                with open ("./profile,rs,interp.txt", "w") as file:
                    print ("# Interpolated profile:", file=file)
                    for q in linspace (qs[0], qs[-1], 3 * len (qs)):
                        print (q, f(p(q)), h(p(q)), e(p(q)), file=file)


exchange (s)


