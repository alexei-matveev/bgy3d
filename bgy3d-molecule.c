/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.15 2007-04-23 16:55:13 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-molecule.h"

#define rab 1.0

static void ComputeDiatomicAB_g(Vec g, Vec g0, Vec dg);

typedef struct BGY3dDiatomicABStruct
{
  DA da, dc;
  Mat fft_mat;

  const ProblemData *PD;

  Vec fa[3],fb[3],fab[3];
  Vec v[3];

  real LJ_paramsa[2], LJ_paramsb[2], LJ_paramsab[2] ; /* sigma and epsilon  */
  real beta, rho;

  real norm_const, c_ab, c_aab;

  Vec ga_ini, gb_ini, gab_ini;

  /* Complex Vecs for parallel FFT: */
  Vec fg2_fft[3], g_fft, gfg2_fft; /* complex */

} *BGY3dDiatomicABData;

static BGY3dDiatomicABData BGY3dDiatomicABData_Pair_malloc(const ProblemData *PD)
{
  BGY3dDiatomicABData BDD;
  real interval[2], h[3], r[3], r_s, beta;
  int i[3], x[3], n[3];
  PetscScalar ***(fa_vec[3]),***(fb_vec[3]),***(fab_vec[3]);
  PetscScalar ***gaini_vec, ***gbini_vec, ***gabini_vec;
  real epsilona, epsilonb, epsilonab;
  real sigmaa, sigmab, sigmaab;

  real eps[]={ 1.0, 1.0}, sig[]={ 0.5, 1.0};
  //real eps[]={ 0.1512, 0.046}, sig[]={ 3.1506, 0.4};


  BDD = (BGY3dDiatomicABData) malloc(sizeof(*BDD));


  BDD->LJ_paramsa[0] = eps[0];
  BDD->LJ_paramsa[1] = sig[0];
  epsilona = BDD->LJ_paramsa[0];
  sigmaa = BDD->LJ_paramsa[1];

  BDD->LJ_paramsb[0] = eps[1];
  BDD->LJ_paramsb[1] = sig[1];
  epsilonb = BDD->LJ_paramsb[0];
  sigmab = BDD->LJ_paramsb[1];

  BDD->LJ_paramsab[0] = sqrt(eps[0]*eps[1]);
  BDD->LJ_paramsab[1] = 0.5*(sig[0]+sig[1]);   
  epsilonab = BDD->LJ_paramsab[0];
  sigmaab = BDD->LJ_paramsab[1];

  BDD->beta = PD->beta;
  BDD->rho  = PD->rho;
  beta = PD->beta;
  BDD->norm_const= PD->h[0]*PD->h[1]*PD->h[2]
    /pow(PD->interval[1]-PD->interval[0],3);

  BDD->PD = PD;

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  FOR_DIM
    h[dim]=PD->h[dim];

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &BDD->fft_mat, &BDD->da, &BDD->dc);
  const DA da = BDD->da;

  /* Create global vectors */
  DACreateGlobalVector(da, &(BDD->ga_ini));
  DACreateGlobalVector(da, &(BDD->gb_ini));
  DACreateGlobalVector(da, &(BDD->gab_ini));
  FOR_DIM
    {
      VecDuplicate(BDD->ga_ini, &(BDD->fa[dim]));
      VecDuplicate(BDD->ga_ini, &(BDD->fb[dim]));
      VecDuplicate(BDD->ga_ini, &(BDD->fab[dim]));
      VecDuplicate(BDD->ga_ini, &(BDD->v[dim]));
    }

  FOR_DIM
    {
      VecSet(BDD->fa[dim],0.0);
      VecSet(BDD->fb[dim],0.0);
      VecSet(BDD->fab[dim],0.0);
    }
  VecSet(BDD->ga_ini, 1.0);
  VecSet(BDD->gb_ini, 1.0);
  VecSet(BDD->gab_ini, 1.0);

  DAVecGetArray(da, BDD->ga_ini, &gaini_vec);
  DAVecGetArray(da, BDD->gb_ini, &gbini_vec);
  DAVecGetArray(da, BDD->gab_ini, &gabini_vec);
  FOR_DIM
    {
      DAVecGetArray(da, BDD->fa[dim], &(fa_vec[dim]));
      DAVecGetArray(da, BDD->fb[dim], &(fb_vec[dim]));
      DAVecGetArray(da, BDD->fab[dim], &(fab_vec[dim]));
    }

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0];


	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );
	  gaini_vec[i[2]][i[1]][i[0]] *=
	    exp(-beta* Lennard_Jones( r_s, epsilona, sigmaa));
	  gbini_vec[i[2]][i[1]][i[0]] *=
	    exp(-beta* Lennard_Jones( r_s, epsilonb, sigmab));
	  gabini_vec[i[2]][i[1]][i[0]] *=
	    exp(-beta* Lennard_Jones( r_s, epsilonab, sigmaab));

	   FOR_DIM
	    {

	      fa_vec[dim][i[2]][i[1]][i[0]] +=
		Lennard_Jones_grad( r_s, r[dim], epsilona, sigmaa);
	      fb_vec[dim][i[2]][i[1]][i[0]] +=
		Lennard_Jones_grad( r_s, r[dim], epsilonb, sigmab);
	      fab_vec[dim][i[2]][i[1]][i[0]] +=
		Lennard_Jones_grad( r_s, r[dim], epsilonab, sigmaab);
	    }
	}

  DAVecRestoreArray(da, BDD->ga_ini, &gaini_vec);
  DAVecRestoreArray(da, BDD->gb_ini, &gbini_vec);
  DAVecRestoreArray(da, BDD->gab_ini, &gabini_vec);
  FOR_DIM
    {
      DAVecRestoreArray(da, BDD->fa[dim], &(fa_vec[dim]));
      DAVecRestoreArray(da, BDD->fb[dim], &(fb_vec[dim]));
      DAVecRestoreArray(da, BDD->fab[dim], &(fab_vec[dim]));
    }

  /* Allocate memory for fft */
  FOR_DIM
    DACreateGlobalVector (BDD->dc, &BDD->fg2_fft[dim]);

  DACreateGlobalVector (BDD->dc, &BDD->g_fft);
  DACreateGlobalVector (BDD->dc, &BDD->gfg2_fft);

  return BDD;
}



static void BGY3dDiatomicABData_free(BGY3dDiatomicABData BDD)
{
  MPI_Barrier( PETSC_COMM_WORLD);

  FOR_DIM
    {
      VecDestroy(BDD->fa[dim]);
      VecDestroy(BDD->fb[dim]);
      VecDestroy(BDD->fab[dim]);
      VecDestroy(BDD->v[dim]);
    }

  VecDestroy(BDD->ga_ini);
  VecDestroy(BDD->gb_ini);
  VecDestroy(BDD->gab_ini);

  /* Free memory for fft */
  FOR_DIM
    VecDestroy (BDD->fg2_fft[dim]);

  VecDestroy (BDD->g_fft);
  VecDestroy (BDD->gfg2_fft);

  DADestroy (BDD->da);
  DADestroy (BDD->dc);
  MatDestroy (BDD->fft_mat);

  free(BDD);
}


static void ComputeDiatomicAB_g(Vec g, Vec g0, Vec dg)
{
  int local_size, i;
  PetscScalar *g_vec, *dg_vec;
  // real g_norm;



  VecGetArray( g, &g_vec);
  VecGetArray( dg, &dg_vec);
  VecGetLocalSize(g, &local_size);

  for(i=0; i<local_size; i++)
    g_vec[i] = exp(-dg_vec[i]);

  VecRestoreArray(g, &g_vec);
  VecRestoreArray(dg, &dg_vec);

  /* g=g0*exp(-dg) */
  VecPointwiseMult(g, g, g0);

  /* apply same normalization for all g functions */
/*   VecSum(g, &g_norm); */
/*   VecScale(g, 1./(BDD->norm_const*g_norm)); */
  //PetscPrintf(PETSC_COMM_WORLD,"Hier:: %e\n",BDD->norm_const*g_norm);

/*   VecView(g0,PETSC_VIEWER_STDERR_WORLD); */
/*   VecView(g,PETSC_VIEWER_STDERR_WORLD); */
}


void VectorMovetoZero(BGY3dDiatomicABData BDD, Vec g)
{
   DA da;
  int x[3], n[2];
  PetscScalar ***g_vec;
  real move;

  da = BDD->da;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if(x[0]==0 && x[1]==0 && x[2]==0 )
    {

      DAVecGetArray(da, g, &g_vec);
      move = g_vec[0][0][0];
      DAVecRestoreArray(da, g, &g_vec);
    }
  MPI_Bcast(&move, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD,"%e", move);

  VecShift(g,-move);

}

void Compute_dg_Pair_inter(BGY3dDiatomicABData BDD,
			   Vec f1[3], real sign1, Vec g1a, Vec g1b,
			   Vec f2[3], real sign2, Vec g2a, Vec g2b,
			   Vec dg, Vec dg_help)
{
  DA da;
  int x[3], n[3], i[3], N[3], ic[3];
  real fac, k_fac, L, k, rho, h, sign;

  const ProblemData *PD = BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  Vec *fg2_fft = BDD->fg2_fft;  /* [3] */

  h=PD->h[0]*PD->h[1]*PD->h[2];
  Vec g_fft = BDD->g_fft;
  Vec dg_fft = BDD->gfg2_fft;
  L = PD->interval[1]-PD->interval[0];
  rho = PD->rho;
  fac = L/(2.*M_PI);  /* BDD->f ist nur grad U, nicht F=-grad U  */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* F1*g1a g1b*/
  /************************************************/

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], g1a, f1[dim]);
      if(sign1 != 1.0)
	VecScale(BDD->v[dim], sign1);
      MatMult (BDD->fft_mat, BDD->v[dim], fg2_fft[dim]);
    }

  /* fft(g-1) */
  MatMult (BDD->fft_mat, g1b, g_fft);

  struct {PetscScalar re, im;} ***g_fft_, ***fg2_fft_[3], ***dg_fft_;

  DAVecGetArray (BDD->dc, g_fft, &g_fft_);
  DAVecGetArray (BDD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecGetArray (BDD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft_[i[2]][i[1]][i[0]].re = 0;
	  dg_fft_[i[2]][i[1]][i[0]].im = 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft_[i[2]][i[1]][i[0]].re = 0;//g_fft[0].re*h;
	      dg_fft_[i[2]][i[1]][i[0]].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = h*h*rho*fac/k;
	      /* phase shift factor for x=x+L/2 */
	      sign = cos(M_PI*ic[0])*cos(M_PI*ic[1])*cos(M_PI*ic[2]);

	      FOR_DIM
		dg_fft_[i[2]][i[1]][i[0]].re += ic[dim] * k_fac * sign *
		(fg2_fft_[dim][i[2]][i[1]][i[0]].re * g_fft_[i[2]][i[1]][i[0]].im +
                 fg2_fft_[dim][i[2]][i[1]][i[0]].im * g_fft_[i[2]][i[1]][i[0]].re) ;

	      FOR_DIM
		dg_fft_[i[2]][i[1]][i[0]].im += ic[dim] * k_fac * sign *
		(- fg2_fft_[dim][i[2]][i[1]][i[0]].re * g_fft_[i[2]][i[1]][i[0]].re
		 + fg2_fft_[dim][i[2]][i[1]][i[0]].im * g_fft_[i[2]][i[1]][i[0]].im);

	    }
	}
  DAVecRestoreArray (BDD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BDD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecRestoreArray (BDD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BDD->fft_mat, dg_fft, dg_help);

  VecScale(dg_help, PD->beta/L/L/L);

  VecCopy(dg_help, dg);

  /************************************************/
  /* F2*g2a g2b*/
  /************************************************/

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], g2a, f2[dim]);
      if(sign2 != 1.0)
	VecScale(BDD->v[dim], sign2);
      MatMult (BDD->fft_mat, BDD->v[dim], fg2_fft[dim]);
    }

  /* fft(g-1) */
  MatMult (BDD->fft_mat, g2b, g_fft);

  DAVecGetArray (BDD->dc, g_fft, &g_fft_);
  DAVecGetArray (BDD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecGetArray (BDD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft_[i[2]][i[1]][i[0]].re = 0;
	  dg_fft_[i[2]][i[1]][i[0]].im = 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft_[i[2]][i[1]][i[0]].re = 0;//g_fft[0].re*h;
	      dg_fft_[i[2]][i[1]][i[0]].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = h*h*rho*fac/k;
	      /* phase shift factor for x=x+L/2 */
	      sign = cos(M_PI*ic[0])*cos(M_PI*ic[1])*cos(M_PI*ic[2]);

	      FOR_DIM
		dg_fft_[i[2]][i[1]][i[0]].re += ic[dim] * k_fac * sign *
		(fg2_fft_[dim][i[2]][i[1]][i[0]].re * g_fft_[i[2]][i[1]][i[0]].im +
                 fg2_fft_[dim][i[2]][i[1]][i[0]].im * g_fft_[i[2]][i[1]][i[0]].re);

	      FOR_DIM
		dg_fft_[i[2]][i[1]][i[0]].im += ic[dim] * k_fac * sign *
		(- fg2_fft_[dim][i[2]][i[1]][i[0]].re * g_fft_[i[2]][i[1]][i[0]].re
		 + fg2_fft_[dim][i[2]][i[1]][i[0]].im * g_fft_[i[2]][i[1]][i[0]].im);

	    }
	}
  DAVecRestoreArray (BDD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BDD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecRestoreArray (BDD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BDD->fft_mat, dg_fft, dg_help);

  VecScale(dg_help, PD->beta/L/L/L);

  VecAXPY(dg,1.0, dg_help);
}


/* Compute intramolecular part */
void Compute_dg_Pair_intra(BGY3dDiatomicABData BDD, Vec f[3], Vec g1, Vec g2,
			   Vec dg, Vec dg_help)
{
  DA da;
  int x[3], n[3], i[3], N[3], ic[3];
  real fac, k_fac, L, k, h, beta; // sign;

  const ProblemData *PD = BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  Vec *fg2_fft = BDD->fg2_fft;  /* [3] */

  h=PD->h[0]*PD->h[1]*PD->h[2];
  Vec g_fft = BDD->g_fft;
  Vec dg_fft = BDD->gfg2_fft;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], g1, f[dim]);
      MatMult (BDD->fft_mat, BDD->v[dim], fg2_fft[dim]);
    }

  /* fft(g) */
  MatMult (BDD->fft_mat, g2, g_fft);

  struct {PetscScalar re, im;} ***g_fft_, ***fg2_fft_[3], ***dg_fft_;

  DAVecGetArray (BDD->dc, g_fft, &g_fft_);
  DAVecGetArray (BDD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecGetArray (BDD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft_[i[2]][i[1]][i[0]].re= 0;
	  dg_fft_[i[2]][i[1]][i[0]].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft_[i[2]][i[1]][i[0]].re = 0;//g_fft[0].re*h;
	      dg_fft_[i[2]][i[1]][i[0]].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = beta*fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;

 	      FOR_DIM
		dg_fft_[i[2]][i[1]][i[0]].re += ic[dim] * k_fac *
		(h * fg2_fft_[dim][i[2]][i[1]][i[0]].im * sin(k)/k);

	      FOR_DIM
		dg_fft_[i[2]][i[1]][i[0]].im += ic[dim] * k_fac *
		(-h * fg2_fft_[dim][i[2]][i[1]][i[0]].re * sin(k)/k);
            }
	}
  DAVecRestoreArray (BDD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BDD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecRestoreArray (BDD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BDD->fft_mat, dg_fft, dg_help);

  VecScale(dg_help, 1./L/L/L);

  VecCopy(dg_help,dg);
}



/* Compute intramolecular part */
static void Compute_dg_Pair_intra_ln(BGY3dDiatomicABData BDD, Vec g, real sign, Vec dg, Vec dg_help)
{
  DA da;
  int x[3], n[3], i[3], N[3], ic[3], local_size;
  real L, k, h;
  PetscScalar *g_vec;


  const ProblemData *PD = BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  Vec g_fft = BDD->g_fft;
  Vec dg_fft = BDD->gfg2_fft;
  L = PD->interval[1]-PD->interval[0];
  /* fac = L/(2.*M_PI); /\* siehe oben ... *\/ */

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* F(g) */
  /************************************************/

  MatMult (BDD->fft_mat, g, g_fft);

  struct {PetscScalar re, im;} ***g_fft_, ***dg_fft_;

  DAVecGetArray (BDD->dc, g_fft, &g_fft_);
  DAVecGetArray (BDD->dc, dg_fft, &dg_fft_);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft_[i[2]][i[1]][i[0]].re= 0;
	  dg_fft_[i[2]][i[1]][i[0]].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft_[i[2]][i[1]][i[0]].re = g_fft_[0][0][0].re * h; /* FIXME? */
	      dg_fft_[i[2]][i[1]][i[0]].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      //k_fac = beta*fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;

	      /* + should be correct ??? */
	      dg_fft_[i[2]][i[1]][i[0]].re += h*g_fft_[i[2]][i[1]][i[0]].re*sin(k)/k;

		}
	}
  DAVecRestoreArray (BDD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BDD->dc, dg_fft, &dg_fft_);

  MatMultTranspose (BDD->fft_mat, dg_fft, dg_help);

  VecScale(dg_help, 1./L/L/L);

  /* ln(g) */
  VecGetArray( dg_help, &g_vec);
  VecGetLocalSize(dg_help, &local_size);

  for (int ijk = 0; ijk < local_size; ijk++)
    g_vec[ijk] = sign * log (g_vec[ijk]);

  VecRestoreArray(dg_help, &g_vec);
  /******************************/

  VecAXPY(dg, 1.0, dg_help);
}




/* Compute normalization condition */
void Compute_dg_Pair_normalization_intra(BGY3dDiatomicABData BDD, Vec g,
					 Vec dg, Vec dg_help)
{
  DA da;
  int x[3], n[3], i[3], N[3], ic[3];
  real L, k, h;

  const ProblemData *PD = BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];


  h=PD->h[0]*PD->h[1]*PD->h[2];
  Vec g_fft = BDD->g_fft;
  Vec dg_fft = BDD->gfg2_fft;
  L = PD->interval[1]-PD->interval[0];
  /* fac = L/(2.*M_PI); /\* siehe oben ... *\/ */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* fft(g/t) */
  MatMult (BDD->fft_mat, g, g_fft);

  struct {PetscScalar re, im;} ***g_fft_, ***dg_fft_;

  DAVecGetArray (BDD->dc, g_fft, &g_fft_);
  DAVecGetArray (BDD->dc, dg_fft, &dg_fft_);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft_[i[2]][i[1]][i[0]].re= 0;
	  dg_fft_[i[2]][i[1]][i[0]].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft_[i[2]][i[1]][i[0]].re = h * g_fft_[0][0][0].re; /* FIXME! */
	      dg_fft_[i[2]][i[1]][i[0]].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      /* k_fac = beta*fac/k; */
	      k = 2.0*M_PI*sqrt(k)*rab/L;

	      /* + oder - hier ??? */
	      dg_fft_[i[2]][i[1]][i[0]].re += h*g_fft_[i[2]][i[1]][i[0]].re*sin(k)/k;
	    }
	}
  DAVecRestoreArray (BDD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BDD->dc, dg_fft, &dg_fft_);

  MatMultTranspose (BDD->fft_mat, dg_fft, dg_help);

  VecScale(dg_help, 1./L/L/L);

  VecCopy(dg_help,dg);
}

/* Compute normalization condition */
void Compute_dg_Pair_normalization(BGY3dDiatomicABData BDD, Vec g1, Vec g2,
				   Vec dg, Vec dg_help)
{
  DA da;
  int x[3], n[3], i[3], N[3], ic[3];
  real L, k, h, sign;

  const ProblemData *PD = BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  Vec fg2_fft = BDD->fg2_fft[0]; /* FIXME? */

  h=PD->h[0]*PD->h[1]*PD->h[2];
  Vec g_fft = BDD->g_fft;
  Vec dg_fft = BDD->gfg2_fft;
  L = PD->interval[1]-PD->interval[0];
  /* fac = L/(2.*M_PI); /\* siehe oben ... *\/ */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* fft(g) */
  MatMult (BDD->fft_mat, g1, g_fft);

  MatMult (BDD->fft_mat, g2, fg2_fft);

  struct {PetscScalar re, im;} ***g_fft_, ***fg2_fft_, ***dg_fft_;

  DAVecGetArray (BDD->dc, g_fft, &g_fft_);
  DAVecGetArray (BDD->dc, dg_fft, &dg_fft_);
  DAVecGetArray (BDD->dc, fg2_fft, &fg2_fft_);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft_[i[2]][i[1]][i[0]].re= 0;
	  dg_fft_[i[2]][i[1]][i[0]].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft_[i[2]][i[1]][i[0]].re = h * h * g_fft_[0][0][0].re * fg2_fft_[0][0][0].re; /* FIXME! */
	      dg_fft_[i[2]][i[1]][i[0]].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      /* k_fac = beta*fac/k; */
	      k = 2.0*M_PI*sqrt(k)*rab/L;
	      /* phase shift factor for x=x+L/2 */
	      sign = cos(M_PI*ic[0])*cos(M_PI*ic[1])*cos(M_PI*ic[2]);

	      dg_fft_[i[2]][i[1]][i[0]].re += h*h*sign*
		(g_fft_[i[2]][i[1]][i[0]].re*fg2_fft_[i[2]][i[1]][i[0]].re -
		 g_fft_[i[2]][i[1]][i[0]].im*fg2_fft_[i[2]][i[1]][i[0]].im );

	      dg_fft_[i[2]][i[1]][i[0]].im += h*h*sign*
		(g_fft_[i[2]][i[1]][i[0]].im*fg2_fft_[i[2]][i[1]][i[0]].re +
		 g_fft_[i[2]][i[1]][i[0]].re*fg2_fft_[i[2]][i[1]][i[0]].im );

	    }
	}
  DAVecRestoreArray (BDD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BDD->dc, dg_fft, &dg_fft_);
  DAVecRestoreArray (BDD->dc, fg2_fft, &fg2_fft_);

  MatMultTranspose (BDD->fft_mat, dg_fft, dg_help);

  VecScale(dg_help, 1./L/L/L/L/L/L);

  VecCopy(dg_help,dg);
}


static void Solve_Normalization_old(BGY3dDiatomicABData BDD, Vec ga, Vec gb, Vec gab, Vec gba,
			 Vec ta, Vec tb, Vec tab, Vec tba, Vec dg, Vec dg_help,
			 real norm_tol)
{
  int i;
  real tab_norm, ta_norm, tb_norm; // a=0.2, g_norm, h;

  VecCopy(ga, ta);
  VecCopy(gb, tb);
  VecCopy(gab, tab);


  Compute_dg_Pair_normalization_intra( BDD, gab, dg, dg_help);
    /* ta */
  VecPointwiseDivide(ta, ga, dg);

  /* tb */
  VecPointwiseDivide(tb, gb, dg);

  /* tab */
  Compute_dg_Pair_normalization_intra( BDD, gb, dg, dg_help);
  VecPointwiseDivide(tab, gab, dg);

  /* tba */
/*   Compute_dg_Pair_normalization_intra( BDD, ga, dg, dg_help); */
/*   VecPointwiseDivide(tba, gba, dg); */

/*   VecView(ta,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

  return ;


  VecSet(ta, 1.0);
  VecSet(tb, 1.0);
  VecSet(tab, 1.0);

  for(i=0; i<200; i++)
    {
      /* tab */
/*       Compute_dg_Pair_normalization_intra( BDD, gb, tb, dg, dg_help); */
/*       VecWAXPY(dg_help, -1.0, dg, tab); */
/*       VecNorm(dg_help, NORM_INFINITY, &tab_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD,"\titer %d: t function norms: tab= %e  ",   */
/* 		  i+1, tab_norm); */
/*       VecAXPBY(tab, a, (1-a), dg); */
      /* ta */
/*       Compute_dg_Pair_normalization_intra( BDD, gab, tab, dg, dg_help); */
/*       VecWAXPY(dg_help, -1.0, dg, ta); */
/*       VecNorm(dg_help, NORM_INFINITY, &ta_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD,"ta= %e  ",   */
/* 		  ta_norm); */
/*       VecAXPBY(ta, a, (1-a), dg); */
      /* tb */
/*       Compute_dg_Pair_normalization_intra( BDD, gab, tab, dg, dg_help); */
/*       VecWAXPY(dg_help, -1.0, dg, tb); */
/*       VecNorm(dg_help, NORM_INFINITY, &tb_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD,"tb= %e  \n",   */
/* 		  tb_norm); */
/*       VecAXPBY(tb, a, (1-a), dg); */
      /* tab */
/*       Compute_dg_Pair_normalization_intra( BDD, ga, ta, dg, dg_help); */
/*       VecWAXPY(dg_help, -1.0, dg, tab); */
/*       VecNorm(dg_help, NORM_INFINITY, &tab_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD,"tab= %e  \n",   */
/* 		  i+1, tab_norm); */
/*       VecAXPBY(tab, a, (1-a), dg); */



      if( tab_norm < norm_tol && ta_norm < norm_tol && tb_norm < norm_tol)
	break;
    }

  VecPointwiseDivide(ta, ga, tb);
  VecPointwiseDivide(tb, ga, tb);
  VecPointwiseDivide(tab, gab, tab);

}


void Compute_dg_Pair_normcorrection(BGY3dDiatomicABData BDD, Vec dg, Vec g)
{
  Vec scratch, re, re_2;

  re = BDD->v[0];
  re_2 = BDD->v[1];
  scratch = BDD->v[2];

  VecSet(re_2, 1.0);

  Compute_dg_Pair_normalization_intra( BDD, g, re, scratch);


/*   VecPointwiseMult(re_2, re, re); */
/*   VecAXPY(re, -1.0, re_2); */


  /* dg/(re+re^2) */
  VecPointwiseDivide(dg, dg, re);
}


/* solve with product ansatz g=g0*dg */
Vec BGY3d_solve_DiatomicAB(const ProblemData *PD, Vec g_ini)
{
  BGY3dDiatomicABData BDD;
  Vec g0a, g0b, g0ab, dga, dgb, dgab, dg_new, dg_new2, f, ga, gb, gab;
  Vec dgba, gba;
  Vec ta, tb, tab, tba;
  int iter;
  PetscScalar dga_norm, dgb_norm, dgab_norm;

  PetscScalar dgba_norm;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving Molecular BGY3d equation with Fourier ansatz...\n");

  BDD = BGY3dDiatomicABData_Pair_malloc(PD);

  /* read BGY3dDiv specific things from command line */
  /* Mixing parameter */
  const real a = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  /*********************************/

  DACreateGlobalVector(BDD->da, &ga);
  DACreateGlobalVector(BDD->da, &gb);
  DACreateGlobalVector(BDD->da, &gab);
  DACreateGlobalVector(BDD->da, &dga);
  DACreateGlobalVector(BDD->da, &dgb);
  DACreateGlobalVector(BDD->da, &dgab);
  DACreateGlobalVector(BDD->da, &dg_new);
  DACreateGlobalVector(BDD->da, &dg_new2);
  DACreateGlobalVector(BDD->da, &f);

  DACreateGlobalVector(BDD->da, &dgba);
  DACreateGlobalVector(BDD->da, &gba);

  DACreateGlobalVector(BDD->da, &ta);
  DACreateGlobalVector(BDD->da, &tb);
  DACreateGlobalVector(BDD->da, &tab);
  DACreateGlobalVector(BDD->da, &tba);

/*   VecDuplicate(BDD->ga_ini, &ga); */
/*   VecDuplicate(BDD->ga_ini, &gb); */
/*   VecDuplicate(BDD->ga_ini, &gab); */
/*   VecDuplicate(BDD->ga_ini, &dga); */
/*   VecDuplicate(BDD->ga_ini, &dgb); */
/*   VecDuplicate(BDD->ga_ini, &dgab); */
/*   VecDuplicate(BDD->ga_ini, &dg_new); */
/*   VecDuplicate(BDD->ga_ini, &f); */

  g0a=BDD->ga_ini;
  g0b=BDD->gb_ini;
  g0ab=BDD->gab_ini;

  /* set initial guess*/
  VecSet(dga,0.0);
  VecSet(dgb,0.0);
  VecSet(dgab,0.0);
  VecSet(dg_new,0.0);

  VecSet(dgba,0.0);

  /* g=g0+exp(-dg) */
  ComputeDiatomicAB_g( ga, g0a, dga);
  ComputeDiatomicAB_g( gb, g0b, dgb);
  ComputeDiatomicAB_g( gab, g0ab, dgab);

  ComputeDiatomicAB_g( gba, g0ab, dgba);

/*   VecScale(gb, BDD->c_ab); */
/*   VecScale(gab, BDD->c_aab); */



  for(iter=0; iter<max_iter; iter++)
    {

      Solve_Normalization_old( BDD,  ga,  gb,  gab, gba, ta, tb, tab, tba, dg_new, f,
			       norm_tol);

      /* f=integral(g) */
      if (1)                    /* kflg was set with -pair */
	{



	   /* g_ab */
	  Compute_dg_Pair_inter(BDD, BDD->fa, 1, ga, gab, BDD->fab, 1, gab, gb,
	  			dg_new, f);


	  Compute_dg_Pair_intra_ln(BDD, tb, -1.0, dg_new, f);

	  Compute_dg_Pair_intra(BDD, BDD->fa, ga, gb, dg_new2, f);
	  Compute_dg_Pair_normcorrection(BDD, dg_new2, ta);
	  VecAXPY(dg_new, 1.0, dg_new2);



	  VecWAXPY(f, -1.0, dg_new, dgab);
	  VecNorm(f, NORM_INFINITY, &dgab_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms: ab= %e  ",
		      iter+1, dgab_norm);
	  VecAXPBY(dgab, a, (1-a), dg_new);
	  ComputeDiatomicAB_g( gab, g0ab, dgab);



	   /* g_ba */
/* 	  Compute_dg_Pair_inter(BDD, BDD->fb, 1, gb, gab, BDD->fab, 1, gab, ga,  */
/* 	  			dg_new, f); */

/* 	  Compute_dg_Pair_intra_ln(BDD, ta, -1.0, dg_new, f); */

/* 	  Compute_dg_Pair_intra(BDD, BDD->fb, gb, ga, dg_new2, f); */
/* 	  Compute_dg_Pair_normcorrection(BDD, dg_new2, tb);  */


/* 	  VecAXPY(dg_new, 1.0, dg_new2); */

/* 	  VecWAXPY(f, -1.0, dg_new, dgab); */
/* 	  VecNorm(f, NORM_INFINITY, &dgab_norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms: ab= %e  ",   */
/* 		      iter+1, dgab_norm); */
/* 	  VecAXPBY(dgab, a, (1-a), dg_new); */
/* 	  ComputeDiatomicAB_g( gab, g0ab, dgab); */

	  /**********************************************/
	  Solve_Normalization_old( BDD,  ga,  gb,  gab, gba, ta, tb, tab, tba, dg_new, f,
			   norm_tol);
	  /**********************************************/

	  /* g_a */
 	  Compute_dg_Pair_inter(BDD, BDD->fa, 1, ga, ga, BDD->fab, 1, gab, gab,
				dg_new, f);

	  Compute_dg_Pair_intra_ln(BDD, tab, -1.0, dg_new, f);

	  Compute_dg_Pair_intra(BDD, BDD->fab, gab, gab, dg_new2, f);
	  Compute_dg_Pair_normcorrection(BDD, dg_new2, tab);
	  VecAXPY(dg_new, 1.0, dg_new2);

	  VecWAXPY(f, -1.0, dg_new, dga);
	  VecNorm(f, NORM_INFINITY, &dga_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"a= %e  ", dga_norm);
	  VecAXPBY(dga, a, (1-a), dg_new);
	  ComputeDiatomicAB_g( ga, g0a, dga);


	  /* g_b */
	  Compute_dg_Pair_inter(BDD, BDD->fb, 1, gb, gb, BDD->fab, 1, gab, gab,
				dg_new, f);

	  Compute_dg_Pair_intra_ln(BDD, tab, -1.0, dg_new, f);

	  Compute_dg_Pair_intra(BDD, BDD->fab, gab, gab, dg_new2, f);
	  Compute_dg_Pair_normcorrection(BDD, dg_new2, tab);
	  VecAXPY(dg_new, 1.0, dg_new2);


	  VecWAXPY(f, -1.0, dg_new, dgb);
	  VecNorm(f, NORM_INFINITY, &dgb_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"b= %e\n", dgb_norm);
	  VecAXPBY(dgb, a, (1-a), dg_new);
	  ComputeDiatomicAB_g( gb, g0b, dgb);





	}
      else {
        // nothing
      }

      /* FIXME: dgba_norm may be used uninitialized here: */
      if(dga_norm<=norm_tol &&  dgb_norm<=norm_tol && dgab_norm<=norm_tol &&
	 dgba_norm<=norm_tol)
	break;
    }

  /* output */
  bgy3d_save_vec ("g00.bin", ga);
  bgy3d_save_vec ("g11.bin", gb);
  bgy3d_save_vec ("g01.bin", gab);

  VecDestroy(ga);
  VecDestroy(gb);
  VecDestroy(gab);
  VecDestroy(dga);
  VecDestroy(dgb);
  VecDestroy(dgab);
  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(f);

  VecDestroy(dgba);
  VecDestroy(gba);
  VecDestroy(ta);
  VecDestroy(tb);
  VecDestroy(tab);
  VecDestroy(tba);

  // ExtractAxis(BDD, g, 0);

  BGY3dDiatomicABData_free(BDD);

  return PETSC_NULL;
}
