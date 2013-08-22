/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3d.h,v 1.66 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdbool.h>            /* bool, true, false */
#include <float.h>              /* DBL_MAX */

#include "petsc.h"
#define VERSION(major, minor) ((major) * 10000 + (minor) * 100)
#define PETSC_VERSION VERSION(PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR)

#if PETSC_VERSION >= VERSION(3, 2)
#  include "petscdmda.h"        /* Vec, Mat, DA, ... */
#else
#  include "petscda.h"          /* Vec, Mat, DA, ... */
#endif

#include "petscdmmg.h"          /* KSP, ... */

/* FIXME: PETSC 3.2 */
#if PETSC_VERSION >= VERSION(3, 2)
typedef DM DA;
#  define VecLoadIntoVector(viewer, vec) VecLoad (vec, viewer)
#  define STENCIL_TYPE          DMDA_STENCIL_STAR
#  define BOUNDARY_TYPE         DMDA_BOUNDARY_PERIODIC
/* before PETSC 3.2 */
#else
#  define STENCIL_TYPE          DA_STENCIL_STAR
#  define BOUNDARY_TYPE         DA_XYZPERIODIC
#  define VecDestroy(x)         (VecDestroy)(*(x))
#  define MatDestroy(x)         (MatDestroy)(*(x))
#  define KSPDestroy(x)         (KSPDestroy)(*(x))
#  define SNESDestroy(x)        (SNESDestroy)(*(x))
#  define PetscViewerDestroy(x) (PetscViewerDestroy)(*(x))
#  define DMDestroy(x)          (DADestroy)(*(x))
#  define VecScatterDestroy(x)  (VecScatterDestroy)(*(x))
#  define ISDestroy(x)          (ISDestroy)(*(x))
#  define PCDestroy(x)          (PCDestroy)(*(x))
#  define DMDAGetCorners        DAGetCorners
#  define DMDACreate3d          DACreate3d
#  define DMCreateGlobalVector  DACreateGlobalVector
#  define DMDAGetInfo           DAGetInfo
#  define DMGetMatrix           DAGetMatrix
#  define DMDAVecGetArray       DAVecGetArray
#  define DMDAVecRestoreArray   DAVecRestoreArray
#  define PetscBool             PetscTruth
#  define PETSC_BOOL            PETSC_TRUTH
#  define PetscOptionsBool      PetscOptionsTruth
/*
  DAGet/RestoreGlobalVector()  is a function-like  macro in  Petsc 3.1
  that expands  to a call to DMGet/RestoreGlobalVector()  with DA cast
  to DM:
*/
#  if PETSC_VERSION < VERSION(3, 1)
#    define DMGetGlobalVector     DAGetGlobalVector
#    define DMRestoreGlobalVector DARestoreGlobalVector
#  else
#    define DMGetGlobalVector(da, x)     (DMGetGlobalVector)((DM)(da), x)
#    define DMRestoreGlobalVector(da, x) (DMRestoreGlobalVector)((DM)(da), x)
#  endif
#endif


/*
 * C99 standards removes M_PI from math.h, the digits were copied from
 * /usr/include/math.h:
 */
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#ifndef BGY3d_H
#define BGY3d_H

/* This macro is used for both reals and integers: */
#define SQR(a)   ((a)*(a))

#define FOR_DIM for(int dim = 0; dim < 3; dim++)

/*
  Distances are  measured in Angstroms  in this code. Site  charges in
  natural units e. The energies appear to be measured in kcal / mol so
  that, for example, the room temperature

    T = 298.15 K = 0.0256 eV = 0.592 kcal [/ mol]

  corresponds to
                           -1
    β = 1 / T = 1.6878 kcal    [* mol]

  Note that  the file  bgy3d-solvents.h quotes 1.6889  instead.  There
  are  several definitions of  calorie in  use.  One  is based  on the
  "International Table calorie" which is 4.1868 J, another is based on
  the  "theormochemical calorie"  which is  4.184 J  exactly  [1].  It
  appears  that at  least  in some  cases  the former,  "International
  Table" definition was assumed in this code.

  [1] http://physics.nist.gov/Pubs/SP811/appenB8.html

  Boltzmann constant assuming the IT-calorie definition:
*/
#define KBOLTZMANN (8.3144621/4186.8) /* kcal/mol/K */

/*
  The interaction energy of two unit charges separated by 1 A is
                   -1
    E = 1 * 1 / 1 A   = 0.529 au = 332 kcal [/ mol]

  The next parameter  appears to have the meaning  of this interaction
  energy of such two unit charges:

    EPSILON0INV = 1 / ε₀

  and is used to  covert electrostatic interaction energies to working
  units.  It  has to be  consistent with other force  field parameters
  defined in bgy3d-solvents.h, notably with Lennard-Jones parameters σ
  and  ε  (FIXME:  so  maybe  it  was  not a  good  idea  to  move  it
  here). These are the original comments:

    You have: e^2/4/pi/epsilon0/angstrom, you want: kcal/avogadro/mol

    => 331.84164

  Again, the code appears to use the IT-calorie to define 1 / ε₀:
*/
#define EPSILON0INV 331.84164 //331.84164

/*
  One Bohr is 0.529 Angstrom.  Bohrs, or rather atomic units, are used
  in QM  codes. This code uses angstroms,  historically.  This literal
  value is used in PG, an alternative is 0.52917725750691647. You will
  NOT need this constant unless  your data comes from external sources
  in atomic units.
*/
#define BOHR 0.52917706

/*
  Pick one harmonic  among the many aliased ones  with the integer FFT
  frequencies k + m * N  and arbitrary m.  Effectively this is used to
  take  non-negative frequencies  k for  0 <=  k <=  N/2  and negative
  frequencies k - N for N/2 < k < N:
*/
#define KFREQ(k, N)  (((k) <= (N) / 2) ? (k) : ((k) - (N)))

/*
  Returns (-1)^k = cos  (πk), that is 1 for even and  -1 for odd k. It
  also needs to work as intended for negative frequencies k:
*/
#define COSSIGN(k)  (((k) % 2) ? -1 : 1)


/* GCC extensions: */
#ifdef __GNUC__
#define GCC_VERSION VERSION(__GNUC__, __GNUC_MINOR__)
#endif

#if GCC_VERSION > VERSION(3, 0)
#define likely(x)    __builtin_expect (!!(x), 1)
#define unlikely(x)  __builtin_expect (!!(x), 0)
#define pure         __attribute__ ((const))
#define deprecated   __attribute__ ((__deprecated__))
#else
#define likely(x)    (x)
#define unlikely(x)  (x)
#define pure                    /* pure */
#define deprecated              /* deprecated */
#endif

#if GCC_VERSION > VERSION(4, 3)
static inline void assert_is_null (void *x)
{
  void **y = (void**) x;
  assert (*y == NULL);
}
#define local __attribute__((cleanup(assert_is_null)))
#else
#define local                   /* local */
#endif

/* Initialized (in the C-speak) in bgy3d.c. See comments there: */
extern int verbosity;

/*
  Much of the code uses "real" to define floating point numbers. There
  is no  guarantee that everything  will work for real  being anything
  else than double.   Was a #define real double  in the original code:
*/
typedef double real;

/* There are a few closures available for HNC-like methods: */
typedef enum {CLOSURE_HNC, CLOSURE_KH, CLOSURE_PY} ClosureEnum;


/* Keep  this  in  sync   with  the  foreign  type  (problem_data)  in
   rism.f90: */
typedef struct ProblemData
{
  /*
    These three are redundant and should fulfil the relation

      L[i] = N[i] * h[i]

    modulo floating point arithmetics, as usual.
  */
  int N[3];                     /* global grid size */
  real L[3];                    /* box size */
  real h[3];                    /* mesh width */
  real beta;                    /* inverse temperature, 1/kT */
  real rho;                     /* solvent density */

  /* Other staff that was retrieved by the solvers themselves from the
     (Petsc) environment: */
  real lambda;                  /* Mixing parameter. */
  real damp;                    /* Scaling factor. */
  int max_iter;                 /* Maximal number of iterations. */
  real norm_tol;                /* Convergence threshold. */
  ClosureEnum closure;          /* HNC, KH, or PY */
} ProblemData;


/* Accessors for often used derived properties: */
static inline real
volume (const ProblemData *PD)
{
  return PD->L[0] * PD->L[1] * PD->L[2];
}


static inline real
volume_element (const ProblemData *PD)
{
  return PD->h[0] * PD->h[1] * PD->h[2];
}


/* Get  problem data  (e.g.  from  command line)  using bgy3d_getopt_*
   interface: */
ProblemData bgy3d_problem_data (void);

real** Load_Molecule (int *N);
void Molecule_free (real **x_M, int N);


typedef struct State
{
  /*
    These are array  descriptors for real and complex  vectors and the
    FFT matrix.  The  data distribution of Petsc vectors  and FFTW MPI
    needs to be consistent, so  that these three should be constructed
    accordingly:
  */
  DA da, dc;
  Mat fft_mat;

  /*
    Distributed  arrays (DA) descriptors  offer an  infrastructure for
    temp Vecs.   You can get them with  DAGetGlobalVector() and return
    with DARestoreGlobalVector().  Do NOT VecDestroy() them! Note that
    the Vecs  are created  as needed, but  up to  ~10 of them  will be
    "cached" and  consume the memory until destroyed  together with DA
    object upon DADestroy().

    Code needs both real and complex work vectors.  Up to four complex
    Vecs   are  used   by  bgy3d_pair()   to  offer   work   space  to
    ComputeFFTfromCoulomb(),      by      Compute_dg_inter(),      and
    Compute_dg_intra(). The last one needs the most.
  */

  /* Immutable command line parameters are stored here: */
  const ProblemData *PD;

#ifdef L_BOUNDARY
  Mat dirichlet_mat;
#endif

#ifdef L_BOUNDARY_MG
  DMMG  *dmmg;
  DA da_dmmg;
#endif

#ifdef WITH_EXTRA_SOLVERS
  /*
    BGY3dM 3-site stuff.  These are vector field quantities indexed by
    two site indices. FIXME: get rid of them, see the 2-site version.

    In much of the  code one refers to the two sites  by literal H and
    O, though the actual sites may  be different.  Let us stick to the
    convention that data  for H is stored in first- and  data for O is
    stored  in  second  position  of  an array.   For  symmetric  pair
    quantities  one  may chose  to  store them  in  a  2x2 array,  say
    g2[2][2] with a constrain that g2[0][1] == g2[1][0].
  */
  Vec fs_g2_fft[2][2][3];       /* complex */
  Vec fl_g2_fft[2][2][3];       /* complex */

  /* Newton stuff */
  Vec wHO_fft, wHH_fft;         /* complex */

  DA da_newton, da_newtonF;
  Vec gH, gHO, gO;
  Vec dgH, dgHO, dgO;
  Vec f, f2, f3, f4;
  Vec pre;
#endif
} State;

State* bgy3d_state_make (const ProblemData *PD);
void bgy3d_state_destroy (State *BHD);
void bgy3d_problem_data_print (const ProblemData *PD);

void bgy3d_comm_allreduce (int n, real x[n]);


/* Sum of an integer array: */
static inline int
isum (int n, const int x[n])
{
  int s = 0;

  for (int i = 0; i < n; i++)
    s += x[i];

  return s;
}


/* Sum of an real array: */
static inline real
sum (int n, const real x[n])
{
  real s = 0.0;

  for (int i = 0; i < n; i++)
    s += x[i];

  return s;
}



/* Returns most  negative number for  zero sized arrays.   Will return
   NaN if there is any in the array. */
static inline double maxval (size_t n, const double x[n])
{
  double max = -DBL_MAX;
  for (size_t i = 0; i < n && !isnan (max); i++)
    {
      /* If x[i] is  NaN comparison will fail and  max will become NaN
         too:  */
      max = x[i] < max ? max : x[i];
      /* We  need  to  break  out  otherwise  in  the  next  iteration
         comparison will  fail again and  max may be overwritten  by a
         valid number ... */
    }
  return max;
}

#endif  /* ifndef BGY3d_H */

