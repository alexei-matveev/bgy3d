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

//#define M_PI 3.141592653589793

#ifndef BGY3d_H
#define BGY3d_H

#define real double
#define SQR(a)   ((a)*(a))
#define FOR_DIM  for(dim=0;dim<3;dim++)

/* which dimension of the vector equation is solved ? */
#define VEC_DIM 1

#define CUTOFF 1.0e+8
#define SHIFT 0.0


#define COSSIGN(i)  (((i)%2)?-1:1)

extern int verbosity;

typedef struct ProblemData
{
  real interval[2];     /* min and max of the domain: 3d-box*/
  real h[3];               /* mesh width */
  real beta;            /* 1/kT */
  real rho;             /* density */
  int N[3], N3;                /* global Grid size */

  real g_xm;            /* g^(N_M)(x_M) */


  /* Parallel stuff */
  int id;                   /* id of this process */
  int np;                   /* number of processes */
  int n[3];                 /* local grid size */
  int nbr_right[3], nbr_left[3];  /* neighboring processes */

} *PData;

/* typedef struct BGY3dField */
/* { */
/*   PetscScalar d[3]; */
/* }Field; */

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
  int vec_dim;        /* Dimension of equation */
  DA da;
  Vec x;             /* grid in real space */
  Vec force;         /* force from all molecule atoms */
  Vec force_single;  /* simple force between 2 atoms */
  Vec Ftimesg2 ;      /* force*g2 */
  FFT_DATA *Ftimesg2_fft;

  Mat M;               /* Matrix for FD-Approximation */
  Vec boundary;        /* Vector for right boundary: g=1 */
  void *LJ_params;   /* sigma and epsilon  */


  Vec v1,v2, v3;        /* Vectors for intermediate results */
  Vec pre;
#ifdef MATPRECOND
  MatPrecond MP;
#endif

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;


  PData PD;            /* pointer to ProblemData */

} *BGY3dParameter;
BGY3dParameter BGY3dParameter_malloc(PData PD, int vdim);
void BGY3dParameter_free(BGY3dParameter params);

typedef struct BGY3dVecStruct
{
  BGY3dParameter params[3];
  Vec fl[3];

} *BGY3dParameterVec;
BGY3dParameterVec BGY3dParameterVec_malloc(PData PD);
void BGY3dParameterVec_free(BGY3dParameterVec par_vec);



typedef struct HNC3dDataStruct
{
  DA da;
  Vec pot;
  Vec h_ini;
  void *LJ_params;   /* sigma and epsilon  */
  real beta, rho;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;

  /* things for arbitrary molecule shape */
  Vec c, v;
  FFT_DATA *c_fft, *h_fft, *ch_fft;



  PData PD;            /* pointer to ProblemData */
} *HNC3dData;

HNC3dData HNC3dData_malloc(PData PD);
void HNC3dData_free(HNC3dData HD);

typedef struct HNC3dField
{
  PetscScalar h, c;
}HNCField;


typedef struct HNC3dNewtonStruct
{
  DA da, da1;
  Vec pot;
  void *LJ_params;   /* sigma and epsilon  */
  real beta, rho;
  Vec pre;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *c_fft, *h_fft, *ch_fft;

  PData PD;            /* pointer to ProblemData */
} *HNC3dNewtonData;

HNC3dNewtonData HNC3dNewtonData_malloc(PData PD);
void HNC3dNewtonData_free(HNC3dNewtonData HD);


/* functions */
real Lennard_Jones(real r, void *LJ_params);
void PData_CreateParallel(PData PD);
real** Load_Molecule(int *N);
void Molecule_free( real **x_M, int N);
void ComputeMatrixStencil(PData PD, DA da, Mat M, int vdim);
Vec BGY3d_solve(PData PD, Vec g_ini, int vec_dim);
void CreateInitialGuess(BGY3dParameter params, Vec g);
void CreateInitialGuessFromg2(BGY3dParameter params, Vec g);
int start_debugger(void );
PetscErrorCode Compute_F(SNES snes, Vec g, Vec f, void *pa);
PetscErrorCode Compute_F_Kirkwood(SNES snes, Vec g, Vec f, void *pa);
PetscErrorCode Compute_J(SNES snes, Vec g, Mat *A, Mat *B, MatStructure *flag,
			 void *pa);
FFT_DATA *ComputeFFTfromVec(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft, int x[3], int n[3], real c);
void ComputeVecfromFFT(DA da, struct fft_plan_3d *fft_plan, Vec g,
		       FFT_DATA *g_fft, int x[3], int n[3], real c);
PetscErrorCode Compute_Preconditioner(void *pa,Vec x,Vec y);
void ConvolutionTest(BGY3dParameter params);

#ifdef MATPRECOND
MatPrecond MatPrecond_malloc(BGY3dParameter params);
void MatPrecond_free(MatPrecond MP);
PetscErrorCode Compute_Preconditioner_Mat(void *pa,Vec x,Vec y);
void TestPreconditioner(MatPrecond MP, Vec x, Vec y);
#endif

/* Vec */
Vec BGY3d_vec_solve(PData PD, Vec g_ini, int vdim);
void CreateInitialGuess_vec(BGY3dParameterVec par_vec, Vec g);
PetscErrorCode ComputeVec_F(SNES snes, Vec g, Vec f, void *pa);

/* HNC3d functions */
void Compute_cgfft(HNC3dData HD, FFT_DATA *c_fft, FFT_DATA *cg_fft, int x[3]
		   ,int  n[3], real h[3]);
void Compute_c_HNC(HNC3dData HD, Vec g, Vec c, int x[3], int n[3]);
Vec HNC3d_Solve(PData PD, Vec g_ini, int vdim);
void SetBoundaryValue(HNC3dNewtonData HD, Vec g, int x[3], int  n[3], real c);
Vec HNC3dNewton_solve(PData PD, Vec g_ini, int vdim);
PetscErrorCode ComputeHNC_F(SNES snes, Vec g, Vec f, void *pa);
PetscErrorCode ComputeHNC_Preconditioner(void *pa,Vec x,Vec y);
void VecOutput_hc(HNC3dNewtonData HD, Vec hc, int horc);
void CreateInitialGuess_HNC(HNC3dNewtonData HD, Vec hc);
PetscErrorCode ComputeHNC2_F(SNES snes, Vec h, Vec f, void *pa);
Vec HNC3dNewton2_solve(PData PD, Vec g_ini, int vdim);
PetscErrorCode ComputeHNC2b_F(SNES snes, Vec h, Vec f, void *pa);
Vec HNC3d_Solve_h(PData PD, Vec g_ini, int vdim);





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

  void *LJ_params;   /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini, g_SA;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;



  PData PD;

}*BGY3dDivData;
BGY3dDivData BGY3dDivData_malloc(PData PD, PetscTruth flg);
BGY3dDivData BGY3dDivData_kirk_malloc(PData PD, PetscTruth flg);
void BGY3dDivData_free(BGY3dDivData BDD);

void AssembleMatrix(BGY3dDivData BDD, DA da, Mat M);
real Lennard_Jones_grad(real r, real xr, void *LJ_params);
Vec BGY3dDiv_solve(PData PD, Vec g_ini, int vdim);
Vec BGY3dDiv_solve2(PData PD, Vec g_ini, int vdim);
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
Vec BGY3dDiv_test(PData PD, Vec g_ini, int vdim);
void InitializeTestData(BGY3dDivData BDD, Vec g, real sigma_g, real sigma_K);
void ComputeRHStest(BGY3dDivData BDD, Vec g, Vec rhs, real sigma_g,
		    real sigma_K);

/* bgy3dfourier */
typedef struct BGY3dFourierStruct
{
  DA da;
  Vec f[3];
  Vec v[3];

  void *LJ_params;   /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;



  PData PD;

}*BGY3dFourierData;
BGY3dFourierData BGY3dFourierData_malloc(PData PD);
BGY3dFourierData BGY3dFourierData_kirk_malloc(PData PD);
void BGY3dFourierData_free(BGY3dFourierData BDD);
void ComputeError(Vec gmax, BGY3dFourierData BDDmax, int Nmax, Vec g, BGY3dFourierData BDD, int N);
void ExtractAxis(BGY3dFourierData BDD, Vec g, int axis);
Vec BGY3dDiv_solve_Fourier(PData PD, Vec g_ini, int vdim);
Vec BGY3dDiv_solve_FourierTest(PData PD, Vec g_ini, int vdim);
Vec BGY3d_Convolution_Test(PData PD, Vec g_ini, int vdim);


/*******************************************/
/* Molecule  */
/*******************************************/


/* diatomic Lennard-Jones */


typedef struct BGY3dDiatomicStruct
{
  DA da;
  Vec f[3];
  Vec v[3];

  void *LJ_params;   /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;



  PData PD;

}*BGY3dDiatomicData;
BGY3dDiatomicData BGY3dDiatomicData_Pair_malloc(PData PD);
void BGY3dDiatomicData_free(BGY3dDiatomicData BDD);
void ComputeDiatomic_g(BGY3dDiatomicData BDD, Vec g, Vec g0, Vec dg);
void Compute_dg_Pair(BGY3dDiatomicData BDD, Vec g, Vec dg);
Vec BGY3d_solve_Diatomic(PData PD, Vec g_ini, int vdim);


typedef struct BGY3dDiatomicABStruct
{
  DA da;
  Vec fa[3],fb[3],fab[3];
  Vec v[3];

  void *LJ_paramsa, *LJ_paramsb,*LJ_paramsab ;   /* sigma and epsilon  */
  real beta, rho;

  real norm_const, c_ab, c_aab;

  Vec ga_ini, gb_ini, gab_ini;


  /* Parallel FFT */
  //struct fft_plan_3d *fft_plan;
  fftw_complex *(fg2_fft[3]), *g_fft, *gfg2_fft, *fft_scratch;

  fftwnd_mpi_plan fft_plan_fw, fft_plan_bw;


  PData PD;

}*BGY3dDiatomicABData;
BGY3dDiatomicABData BGY3dDiatomicABData_Pair_malloc(PData PD);
void BGY3dDiatomicABData_free(BGY3dDiatomicABData BDD);
void Compute_dg_Pair_inter(BGY3dDiatomicABData BDD, Vec f1[3], real sign1,
			   Vec g1a, Vec g1b,
			   Vec f2[3], real sign2,
			   Vec g21, Vec g2b, Vec dg, Vec g_help);
void Compute_dg_Pair_intra(BGY3dDiatomicABData BDD, Vec f[3], Vec g1, Vec g2,
			   Vec dg, Vec dg_help);
void Compute_dg_Pair_normalization(BGY3dDiatomicABData BDD, Vec g1, Vec g2,
				   Vec dg, Vec dg_help);
Vec BGY3d_solve_DiatomicAB(PData PD, Vec g_ini, int vdim);
void ComputeDiatomicAB_g(Vec g, Vec g0, Vec dg);
fftw_complex *ComputeFFTfromVec_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
				fftw_complex *g_fft, fftw_complex *work,
				int x[3], int n[3], real c);
void ComputeVecfromFFT_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
			    fftw_complex *g_fft, fftw_complex *work,
			    int x[3], int n[3], real c);



/*******************************************/
/* Water  */
/*******************************************/
typedef struct BGY3dH2OStruct
{
  DA da;
  Vec fH[3],fO[3],fHO[3];
  Vec fH_l[3], fO_l[3], fHO_l[3];
  Vec ucH, ucHO, ucO;
  Vec v[3];
  Vec g2H, g2O, g2HO;

  Vec cH, cHO, cO;


  void *LJ_paramsH, *LJ_paramsO,*LJ_paramsHO ;   /* sigma and epsilon  */
  real beta, rho;
  real rho_H, rho_O;

  Vec gH_ini, gO_ini, gHO_ini;
  real ucH_0, ucO_0, ucHO_0;

  /* Parallel FFT */
  //struct fft_plan_3d *fft_plan;
  fftw_complex *(fg2_fft[3]), *g_fft, *gfg2_fft, *fft_scratch;
  fftw_complex *ucH_fft, *ucO_fft, *ucHO_fft;
  fftw_complex *wHO_fft, *wHH_fft;

  fftwnd_mpi_plan fft_plan_fw, fft_plan_bw;
  int p_id, p_index;

  PData PD;

  /* BGY3dM stuff */
  fftw_complex *(fg2OO_fft[3]), *(fg2HH_fft[3]), *(fg2HO_fft[3]);
  fftw_complex *(fg2OOl_fft[3]), *(fg2HHl_fft[3]), *(fg2HOl_fft[3]);
  fftw_complex *(fO_fft[3]), *(fH_fft[3]);

#ifdef L_BOUNDARY
  Mat M;
  KSP ksp;
  Vec xH, xHO, xO;   /* solutions of the Laplace */
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

}*BGY3dH2OData;

typedef struct H2OdgStruct
{
  PetscScalar dgH, dgO, dgHO;
} H2Odg;

typedef struct H2OSdgStruct
{
  PetscScalar dgH, dgO;
} H2OSdg;

typedef struct H2OSdgFStruct
{
  PetscScalar dgHre, dgHim, dgOre, dgOim;

} H2OSdgF;

void Smooth_Function(BGY3dH2OData BHD, Vec g, real RL, real RR, real shift);
void Zeropad_Function(BGY3dH2OData BHD, Vec g, real ZP, real shift);
BGY3dH2OData BGY3dH2OData_Pair_malloc(PData PD);
void BGY3dH2OData_free(BGY3dH2OData BHD);
Vec BGY3d_solve_2site(PData PD, Vec g_ini, int vdim);
Vec BGY3d_solve_3site(PData PD, Vec g_ini, int vdim);
Vec BGY3d_solve_4site(PData PD, Vec g_ini, int vdim);
// Replaced void pointer as real number in Coulomb_short() and Coulomb_short_grad()
// real Coulomb_short( real r, void *params);
real Coulomb_short( real r, real SQRq);
// real Coulomb_short_grad( real r, real rx, void *params);
real Coulomb_short_grad( real r, real rx, real SQRq);
real Coulomb_long( real r, void *params);
real Coulomb_long_grad( real r, real rx, void *params);
real Coulomb( real r, void *params);
real Coulomb_grad( real r, real rx, void *params);
void ComputeFFTfromCoulomb(BGY3dH2OData BHD, Vec uc, Vec f_l[3], fftw_complex *fft_data, void *LJ_params, real damp);
void ComputeFFTfromCoulombII(BGY3dH2OData BHD, Vec f[3], Vec f_l[3], fftw_complex *fft_data, void *LJ_params, real damp);
void ComputeFFTSoluteII(BGY3dH2OData BHD, Vec ucl , Vec ucs, void *LJ_params,
			real damp, real zpad);
void ComputeH2O_g(Vec g, Vec g0, Vec dg);
void Compute_dg_H2O_inter(BGY3dH2OData BHD,
			  Vec f1[3], Vec f1_l[3], Vec g1a, Vec g1b,
			  fftw_complex *coul1_fft, real rho1, real shift1,
			  Vec f2[3], Vec f2_l[3], Vec g2a, Vec g2b,
			  fftw_complex *coul2_fft, real rho2, real shift2,
			  Vec dg, Vec dg_help);
void Compute_dg_H2O_intra(BGY3dH2OData BHD, Vec f[3], Vec f_l[3], Vec g1, Vec g2,
			  fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help);
void Compute_dg_H2O_intra_ln(BGY3dH2OData BHD, Vec g, real rab, Vec dg, Vec dg_help);
void Compute_dg_H2O_intra_lnII(BGY3dH2OData BHD, Vec g, Vec t, real rab, Vec dg, Vec dg_help);
void Compute_dg_H2O_intra_lnIII(BGY3dH2OData BHD, Vec g, Vec t, real rab, Vec dg, Vec dg_help);
void Solve_NormalizationH2O_small(BGY3dH2OData BHD, Vec gc, real rc, Vec g, Vec t,
				  Vec dg, Vec dg_help, real zpad);
void Solve_NormalizationH2O_smallII(BGY3dH2OData BHD, Vec gc, real rc, Vec g, Vec t,
				  Vec dg, Vec dg_help, real zpad);
Vec BGY3d_SolveNewton_H2O(PData PD, Vec g_ini, int vdim);
Vec BGY3d_SolveNewton_H2OS(PData PD, Vec g_ini, int vdim);
void RecomputeInitialData(BGY3dH2OData BHD, real damp, real damp_LJ);
void VecSetRandom_H2O(Vec g, real mag);
void Compute_dg_H2O_intraII(BGY3dH2OData BHD, Vec f[3], Vec f_l[3], Vec g1, Vec tg,
			    fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help);
void Compute_dg_H2O_intraIII(BGY3dH2OData BHD, Vec f[3], Vec f_l[3], Vec g1, Vec tg,
                             fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help);
void Compute_dg_H2O_normalization_intra(BGY3dH2OData BHD, Vec g, real rab,
					Vec dg, Vec dg_help);
void ImposeBoundaryCondition_Initialize( BGY3dH2OData BHD, real zpad);
void ImposeBoundaryCondition( BGY3dH2OData BHD, Vec g);
real ImposeBoundaryConditionII( BGY3dH2OData BHD, Vec g, real zpad);
void WriteH2ONewtonPlain(BGY3dH2OData BHD, Vec u);
void WriteH2ONewtonSolution(BGY3dH2OData BHD, Vec u);
/* BGY3dM */
void ReadPairDistribution(BGY3dH2OData BHD, char *filename, Vec g2);
Vec BGY3dM_solve_H2O_3site(PData PD, Vec g_ini, int vdim);
Vec BGY3dM_solve_H2O_2site(PData PD, Vec g_ini, int vdim);
void RecomputeInitialFFTs(BGY3dH2OData BHD, real damp, real damp_LJ);
void RecomputeInitialSoluteData(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
real ComputeCharge(BGY3dH2OData BHD, Vec g1, Vec g2);
void Compute_H2O_interS(BGY3dH2OData BHD,
			fftw_complex *(fg2_fft[3]), Vec g, fftw_complex *coul_fft,
			fftw_complex *(fs_fft[3]), real con, real rho, Vec dg_help);
void WriteH2OSNewtonSolution(BGY3dH2OData BHD, Vec u);
void WriteH2OSNewtonPlain(BGY3dH2OData BHD, Vec u);
void EnforceNormalizationCondition(BGY3dH2OData BHD, Vec dgO, Vec dgH, Vec gO, Vec gH);
Vec BGY3d_SolveNewton_H2OSF(PData PD, Vec g_ini, int vdim);
#ifdef L_BOUNDARY
void InitializeLaplaceMatrix(BGY3dH2OData BHD, real zpad);
void InitializeKSPSolver(BGY3dH2OData BHD);
real  ImposeLaplaceBoundary(BGY3dH2OData BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif
#ifdef L_BOUNDARY_MG
void InitializeDMMGSolver(BGY3dH2OData BHD);
real  ImposeLaplaceBoundary(BGY3dH2OData BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif
/* Solute Functions */
void RecomputeInitialSoluteData_Methanol(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
void RecomputeInitialSoluteData_Water(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
void RecomputeInitialSoluteData_CS2(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
void RecomputeInitialSoluteData_ButanoicAcid(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
void RecomputeInitialSoluteData_Hexane(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
void RecomputeInitialSoluteData_HCl(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
#endif

