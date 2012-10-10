/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
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
#include <fftw_mpi.h>

#include "petscsnes.h"
#include "petscda.h"
#include "fft_3d.h"
#include "petscdmmg.h"

/*
 * C99 standards removes M_PI from math.h, the digits were copied from
 * /usr/include/math.h:
 */
#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#ifndef BGY3d_H
#define BGY3d_H

#define real double
#define SQR(a)   ((a)*(a))
#define FOR_DIM for(int dim = 0; dim < 3; dim++)

/*  One Bohr  is 0.53  Angstrom.  Bohrs, or  rather atomic  units, are
    used in QM codes. This code uses historically angstroms: */
#define BOHR 0.52917706 /* <- used in PG, 0.52917725750691647 */
#define CUTOFF 1.0e+8
#define SHIFT 0.0


#define COSSIGN(i)  (((i)%2)?-1:1)

extern int verbosity;

typedef struct ProblemData
{
  real interval[2];             /* min and max of the domain: 3d-box*/
  real h[3];                    /* mesh width */
  real beta;                    /* 1/kT */
  real rho;                     /* density */
  int N[3], N3;                 /* global Grid size */

  /* Parallel stuff */
  int id;                       /* id of this process */
  int np;                       /* number of processes */

  /* Other staff that was retrieved by the solvers themselves from the
     (Petsc) environment: */
  real lambda;                  /* Mixing parameter. */
  real damp;                    /* Scaling factor. */
  int max_iter;                 /* Maximal number of iterations. */
  real norm_tol;                /* Convergence threshold. */
  real zpad;                    /* FIXME: ??? */
} ProblemData;

/* Get  problem data  (e.g.  from  command line)  using bgy3d_getopt_*
   interface: */
ProblemData bgy3d_problem_data (void);

/* functions */
real Lennard_Jones(real r, real epsilon, real sigma);
real Lennard_Jones_grad(real r, real xr, real epsilon, real sigma);

real** Load_Molecule(int *N);
void Molecule_free( real **x_M, int N);

/*******************************************/
/* Water  */
/*******************************************/
typedef struct State
{
  /*
   * In much of the code one refers  to the two sites by literal H and
   * O, though the actual sites may be different.  Let us stick to the
   * convention that data for H is  stored in first- and data for O is
   * stored  in  second  position  of  an array.  For  symmetric  pair
   * quantities  one may  chose  to store  them  in a  2x2 array,  say
   * g2[2][2] with a constrain that g2[0][1] == g2[1][0].
   */
  DA da;
  Vec F[2][2][3];               /* sort range pair force */
  Vec F_l[2][2][3];             /* long range pair force, redundant */
  Vec v[3];                     /* work vectors */
  Vec g2[2][2];                 /* Site-site  distributions. Used only
                                   in the solute-solvent solvers. */

  /* Long-range Coulomb interaction for solvent site pairs. So far the
     pairs differ only by a factor q[i] * q[j]. Maybe we should rather
     store just one? */
  Vec u2[2][2];

  fftw_complex *u2_fft[2][2];   /* The  fourier transform of  u2.  The
                                   same redundancy. */

  Vec cH, cHO, cO;

  /* Pair  interaction parameters.  Not used  with  impurities.  Array
     entries are: sigma, epsilon and charge product. */
  real LJ_paramsH[3], LJ_paramsO[3], LJ_paramsHO[3];

  real beta, rho;          /* inverse temperature, solvent density. */

  /*
   * The solute  field for  each of the  two solvent sites  (scaled by
   * inverse  temperature beta)  is initially  put into  the following
   * array:
   */
  Vec g_ini[2];                 /* Short-range  force   field  of  the
                                   solute   for  H   and  O   in  that
                                   order. */

  real rhos[2];                 /* Site specific density.  Computed as
                                   a solvent  density rho times number
                                   of   sites  of   that  type   in  a
                                   solvent. */

  Vec gHO_ini;                  /* used for pure solvent only */

  /*
    Parallel FFT.

    g_fft appears to be used as a temporary in
    ComputeSoluteDatafromCoulomb*() group of functions.
   */
  fftw_complex *fg2_fft[3], *g_fft, *gfg2_fft, *fft_scratch;
  fftw_complex *wHO_fft, *wHH_fft; /* used for pure solvent only */

  fftwnd_mpi_plan fft_plan_fw, fft_plan_bw;
  int p_id, p_index;

  const ProblemData *PD;

  /* BGY3dM stuff.   These are vector field quantities  indexed by two
     site indices */
  fftw_complex *f_g2_fft[2][2][3];
  fftw_complex *fl_g2_fft[2][2][3];
  fftw_complex *fO_fft[3], *fH_fft[3];

#ifdef L_BOUNDARY
  Mat M;
  KSP ksp;
  Vec x_lapl[2], xHO;           /* solutions of the Laplace */
#endif

#ifdef L_BOUNDARY_MG
  DMMG  *dmmg;
  DA da_dmmg;
#endif


  /* Newton stuff */
  DA da_newton, da_newtonF;
  Vec gH, gHO, gO;
  Vec dgH, dgHO, dgO;
  Vec f, f2, f3, f4;
  Vec pre;
  real zpad;

} State;

#endif  /* ifndef BGY3d_H */

