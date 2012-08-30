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

/* which dimension of the vector equation is solved ? */
/* UNUSED: #define VEC_DIM 1 */

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

  real g_xm;                    /* g^(N_M)(x_M) */

  /* Parallel stuff */
  int id;                        /* id of this process */
  int np;                        /* number of processes */
  int n[3];                      /* local grid size */
  int nbr_right[3], nbr_left[3]; /* neighboring processes */
} ProblemData;

/*
 * FIXME: "Never _ever_ make the  "pointerness" part of the type", (by
 * Linus Torvalds).  Consider converting "PType x" to "Type *x".
 */

#ifdef MATPRECOND
typedef struct MatPrecondStruct
{
  Mat P;
  KSP ksp;
  PC  pc;
} *MatPrecond;
#endif

typedef struct BGY3dParameter
{
  int vec_dim;                  /* Dimension of equation */
  DA da;
  Vec x;                        /* grid in real space */
  Vec force;                    /* force from all molecule atoms */
  Vec force_single;             /* simple force between 2 atoms */
  Vec Ftimesg2 ;                /* force*g2 */
  FFT_DATA *Ftimesg2_fft;

  Mat M;                        /* Matrix for FD-Approximation */
  Vec boundary;                 /* Vector for right boundary: g=1 */
  real LJ_params[2]; /* sigma and epsilon, seems we don't need charge here */


  Vec v1,v2, v3;                /* Vectors for intermediate results */
  Vec pre;
#ifdef MATPRECOND
  MatPrecond MP;
#endif

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;

  ProblemData *PD;
} BGY3dParameter;

BGY3dParameter *BGY3dParameter_malloc(ProblemData *PD, int vdim);
void BGY3dParameter_free(BGY3dParameter *params);

/* functions */
real Lennard_Jones(real r, real epsilon, real sigma);
real Lennard_Jones_grad(real r, real xr, real epsilon, real sigma);

real** Load_Molecule(int *N);
void Molecule_free( real **x_M, int N);
void CreateInitialGuess(BGY3dParameter *params, Vec g);
PetscErrorCode Compute_F(SNES snes, Vec g, Vec f, void *pa);

PetscErrorCode Compute_Preconditioner(void *pa,Vec x,Vec y);

#ifdef MATPRECOND
MatPrecond MatPrecond_malloc(BGY3dParameter *params);
void MatPrecond_free(MatPrecond MP);
PetscErrorCode Compute_Preconditioner_Mat(void *pa,Vec x,Vec y);
void TestPreconditioner(MatPrecond MP, Vec x, Vec y);
#endif

/*******************************************/
/* Water  */
/*******************************************/
typedef struct State
{
  DA da;
  Vec fH[3],fO[3],fHO[3];
  Vec fH_l[3], fO_l[3], fHO_l[3];
  Vec v[3];
  Vec g2H, g2O, g2HO;

  Vec cH, cHO, cO;

  real LJ_paramsH[3], LJ_paramsO[3], LJ_paramsHO[3] ; /* sigma, epsilon and charge(product)  */
  real beta, rho;

  /*
   * The solute  field for  each of the  two solvent sites  (scaled by
   * inverse  temperature beta)  is initially  put into  the following
   * array.
   *
   * In much of the code one refers  to the two sites by literal H and
   * O, though the actual sites may be different.  Let us stick to the
   * convention that data for H is  stored in first- and data for O is
   * stored  in second  position of  an  array. Here  in g_ini[0]  and
   * g_ini[1] or uc[0] and uc[1], respectively.
   */
  Vec g_ini[2];                 /* Short-range  force   field  of  the
                                   solute   for  H   and  O   in  that
                                   order. */

  Vec uc[2];                    /* Long-range Coulomb  field for H and
                                   O in that order. */

  real rhos[2];                 /* Site specific density.  Computed as
                                   a solvent  density rho times number
                                   of   sites  of   that  type   in  a
                                   solvent. */

  Vec gHO_ini;
  Vec ucHO;

  real ucH_0, ucO_0, ucHO_0;

  /*
    Parallel FFT.

    g_fft appears to be used as a temporary in
    ComputeSoluteDatafromCoulomb*() group of functions.
   */
  //struct fft_plan_3d *fft_plan;
  fftw_complex *(fg2_fft[3]), *g_fft, *gfg2_fft, *fft_scratch;
  fftw_complex *ucH_fft, *ucO_fft, *ucHO_fft;
  fftw_complex *wHO_fft, *wHH_fft;

  fftwnd_mpi_plan fft_plan_fw, fft_plan_bw;
  int p_id, p_index;

  ProblemData *PD;

  /* BGY3dM stuff */
  fftw_complex *(fg2OO_fft[3]), *(fg2HH_fft[3]), *(fg2HO_fft[3]);
  fftw_complex *(fg2OOl_fft[3]), *(fg2HHl_fft[3]), *(fg2HOl_fft[3]);
  fftw_complex *(fO_fft[3]), *(fH_fft[3]);

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

