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

typedef struct BGY3dParameterStruct
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
} *BGY3dParameter;
BGY3dParameter BGY3dParameter_malloc(ProblemData *PD, int vdim);
void BGY3dParameter_free(BGY3dParameter params);

typedef struct BGY3dVecStruct
{
  BGY3dParameter params[3];
  Vec fl[3];

} *BGY3dParameterVec;
BGY3dParameterVec BGY3dParameterVec_malloc(ProblemData *PD);
void BGY3dParameterVec_free(BGY3dParameterVec par_vec);



typedef struct HNC3dDataStruct
{
  DA da;
  Vec pot;
  Vec h_ini;
  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;

  /* things for arbitrary molecule shape */
  Vec c, v;
  FFT_DATA *c_fft, *h_fft, *ch_fft;

  ProblemData *PD;
} *HNC3dData;

HNC3dData HNC3dData_malloc(ProblemData *PD);
void HNC3dData_free(HNC3dData HD);

typedef struct HNCField
{
  PetscScalar h, c;
} HNCField;


typedef struct HNC3dNewtonStruct
{
  DA da, da1;
  Vec pot;
  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;
  Vec pre;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *c_fft, *h_fft, *ch_fft;

  ProblemData *PD;
} *HNC3dNewtonData;

HNC3dNewtonData HNC3dNewtonData_malloc(ProblemData *PD);
void HNC3dNewtonData_free(HNC3dNewtonData HD);


/* functions */
real Lennard_Jones(real r, real epsilon, real sigma);
void PData_CreateParallel(ProblemData *PD);
real** Load_Molecule(int *N);
void Molecule_free( real **x_M, int N);
void ComputeMatrixStencil(ProblemData *PD, DA da, Mat M, int vdim);
Vec BGY3d_solve(ProblemData *PD, Vec g_ini, int vec_dim);
void CreateInitialGuess(BGY3dParameter params, Vec g);
void CreateInitialGuessFromg2(BGY3dParameter params, Vec g);
int start_debugger(void );
PetscErrorCode Compute_F(SNES snes, Vec g, Vec f, void *pa);
PetscErrorCode Compute_F_Kirkwood(SNES snes, Vec g, Vec f, void *pa);
PetscErrorCode Compute_J(SNES snes, Vec g, Mat *A, Mat *B, MatStructure *flag,
			 void *pa);

PetscErrorCode Compute_Preconditioner(void *pa,Vec x,Vec y);
void ConvolutionTest(BGY3dParameter params);

#ifdef MATPRECOND
MatPrecond MatPrecond_malloc(BGY3dParameter params);
void MatPrecond_free(MatPrecond MP);
PetscErrorCode Compute_Preconditioner_Mat(void *pa,Vec x,Vec y);
void TestPreconditioner(MatPrecond MP, Vec x, Vec y);
#endif

/* Vec */
Vec BGY3d_vec_solve(ProblemData *PD, Vec g_ini, int vdim);
void CreateInitialGuess_vec(BGY3dParameterVec par_vec, Vec g);
PetscErrorCode ComputeVec_F(SNES snes, Vec g, Vec f, void *pa);

/* HNC3d functions */
void Compute_cgfft(HNC3dData HD, FFT_DATA *c_fft, FFT_DATA *cg_fft, int x[3]
		   ,int  n[3], real h[3]);
void Compute_c_HNC(HNC3dData HD, Vec g, Vec c, int x[3], int n[3]);
Vec HNC3d_Solve(ProblemData *PD, Vec g_ini, int vdim);
void SetBoundaryValue(HNC3dNewtonData HD, Vec g, int x[3], int  n[3], real c);
Vec HNC3dNewton_solve(ProblemData *PD, Vec g_ini, int vdim);
PetscErrorCode ComputeHNC_F(SNES snes, Vec g, Vec f, void *pa);
PetscErrorCode ComputeHNC_Preconditioner(void *pa,Vec x,Vec y);
void VecOutput_hc(HNC3dNewtonData HD, Vec hc, int horc);
void CreateInitialGuess_HNC(HNC3dNewtonData HD, Vec hc);
PetscErrorCode ComputeHNC2_F(SNES snes, Vec h, Vec f, void *pa);
Vec HNC3dNewton2_solve(ProblemData *PD, Vec g_ini, int vdim);
PetscErrorCode ComputeHNC2b_F(SNES snes, Vec h, Vec f, void *pa);
Vec HNC3d_Solve_h(ProblemData *PD, Vec g_ini, int vdim);





/*******************************************/

typedef struct BGY3dDivStruct
{
  DA da;
  Vec f[3];
  Vec ddU;
  Vec boundary;
  Vec v[3], i[3], v2[3];
  Mat M;
  Mat FD[3];

  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini, g_SA;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;

  ProblemData *PD;
} *BGY3dDivData;
BGY3dDivData BGY3dDivData_malloc(ProblemData *PD, PetscTruth flg);
BGY3dDivData BGY3dDivData_kirk_malloc(ProblemData *PD, PetscTruth flg);
void BGY3dDivData_free(BGY3dDivData BDD);

void AssembleMatrix(BGY3dDivData BDD, DA da, Mat M);
real Lennard_Jones_grad(real r, real xr, real epsilon, real sigma);
Vec BGY3dDiv_solve(ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3dDiv_solve2(ProblemData *PD, Vec g_ini, int vdim);
void AssembleFDMatrix(BGY3dDivData BDD, DA da, Mat M, int vdim);
void ComputeIntegralPart(BGY3dDivData BDD, Vec g, Vec f);
void ComputeIntegralPart_kirk(BGY3dDivData BDD, Vec g, Vec f);
void AssembleSystemMatrix(BGY3dDivData BDD, Mat SM, Vec f);
void AssembleSystemMatrix_part2(BGY3dDivData BDD, Mat SM);
void ComputeRHS(BGY3dDivData BDD, Vec b, Vec g0, Vec f);
void ComputeBGY3dDiv_F(BGY3dDivData BDD, Mat SM, Vec g0, Vec dg, Vec g,
		       Vec b, Vec f);
void ComputeRHS2(BGY3dDivData BDD, Vec b);
void ShiftVec(DA da, Vec g, Vec scratch, int N[3]);
void AssembleSystemMatrix_part2b(BGY3dDivData BDD, Mat SM);

/* bgy3dtest */
Vec BGY3dDiv_test(ProblemData *PD, Vec g_ini, int vdim);
void InitializeTestData(BGY3dDivData BDD, Vec g, real sigma_g, real sigma_K);
void ComputeRHStest(BGY3dDivData BDD, Vec g, Vec rhs, real sigma_g,
		    real sigma_K);

/* bgy3dfourier */
typedef struct BGY3dFourierStruct
{
  DA da;
  Vec f[3];
  Vec v[3];

  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;

  ProblemData *PD;
} *BGY3dFourierData;
BGY3dFourierData BGY3dFourierData_malloc(ProblemData *PD);
BGY3dFourierData BGY3dFourierData_kirk_malloc(ProblemData *PD);
void BGY3dFourierData_free(BGY3dFourierData BDD);
void ComputeError(Vec gmax, BGY3dFourierData BDDmax, int Nmax, Vec g, BGY3dFourierData BDD, int N);
void ExtractAxis(BGY3dFourierData BDD, Vec g, int axis);
Vec BGY3dDiv_solve_Fourier(ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3dDiv_solve_FourierTest(ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3d_Convolution_Test(ProblemData *PD, Vec g_ini, int vdim);


/*******************************************/
/* Molecule  */
/*******************************************/


/* diatomic Lennard-Jones */


typedef struct BGY3dDiatomicStruct
{
  DA da;
  Vec f[3];
  Vec v[3];

  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;

  ProblemData *PD;
} *BGY3dDiatomicData;
BGY3dDiatomicData BGY3dDiatomicData_Pair_malloc(ProblemData *PD);
void BGY3dDiatomicData_free(BGY3dDiatomicData BDD);
void ComputeDiatomic_g(BGY3dDiatomicData BDD, Vec g, Vec g0, Vec dg);
void Compute_dg_Pair(BGY3dDiatomicData BDD, Vec g, Vec dg);
Vec BGY3d_solve_Diatomic(ProblemData *PD, Vec g_ini, int vdim);


typedef struct BGY3dDiatomicABStruct
{
  DA da;
  Vec fa[3],fb[3],fab[3];
  Vec v[3];

  real LJ_paramsa[2], LJ_paramsb[2], LJ_paramsab[2] ; /* sigma and epsilon  */
  real beta, rho;

  real norm_const, c_ab, c_aab;

  Vec ga_ini, gb_ini, gab_ini;


  /* Parallel FFT */
  //struct fft_plan_3d *fft_plan;
  fftw_complex *(fg2_fft[3]), *g_fft, *gfg2_fft, *fft_scratch;

  fftwnd_mpi_plan fft_plan_fw, fft_plan_bw;

  ProblemData *PD;
} *BGY3dDiatomicABData;

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

typedef struct H2Odg
{
  PetscScalar dgH, dgO, dgHO;
} H2Odg;

typedef struct H2OSdg
{
  PetscScalar dgH, dgO;
} H2OSdg;

typedef struct H2OSdgF
{
  PetscScalar dgHre, dgHim, dgOre, dgOim;

} H2OSdgF;

Vec BGY3d_solve_4site(ProblemData *PD, Vec g_ini, int vdim);

#ifdef L_BOUNDARY_MG
void InitializeDMMGSolver (State *BHD);
real ImposeLaplaceBoundary (State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif
#endif

