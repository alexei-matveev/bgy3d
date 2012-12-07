/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O.c,v 1.42 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-solvents.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-getopt.h"
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-impure.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-poisson.h"      /* laplace staff */
#include "bgy3d-pure.h"
#include <complex.h>            /* after fftw.h */

static real NORM_REG=1.0e-1;
static real NORM_REG2=1.0e-2;

static void normalization_intra (const State *BHD,
                                 Vec g_fft, /* complex, in */
                                 real rab,
                                 Vec work, /* complex, work */
                                 Vec dg);  /* real, out */

static State *BGY3dH2OData_Pair_malloc (const ProblemData *PD)
{
  State *BHD;

  BHD = (State*) malloc(sizeof(*BHD));

  /****************************************************/
  /* set Lennard-Jones and Coulomb parameters */
  /****************************************************/

  /* water hydrogen */
  BHD->LJ_paramsH[0] = eH;  /* epsilon  */
  BHD->LJ_paramsH[1] = sH;  /* sigma    */
  BHD->LJ_paramsH[2] = SQR(qH); /* charge product */

  /* water oxygen */
  BHD->LJ_paramsO[0] = eO;  /* epsilon  */
  BHD->LJ_paramsO[1] = sO;  /* sigma    */
  BHD->LJ_paramsO[2] = SQR(qO); /* charge product */

  /* water O-H mixed parameters */
  BHD->LJ_paramsHO[0] = sqrt(eH*eO);  /* epsilon  */
  BHD->LJ_paramsHO[1] = 0.5*(sH+sO);  /* sigma    */
  BHD->LJ_paramsHO[2] = qH*qO; /* charge product */

  /****************************************************/


  BHD->PD = PD;

  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  PetscPrintf(PETSC_COMM_WORLD, "Regularization of normalization: NORM_REG= %e\n", NORM_REG);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", PD->beta);
  /******************************/
  PetscPrintf(PETSC_COMM_WORLD,"--------------------------------\n");
  PetscPrintf(PETSC_COMM_WORLD,"eps1= %f \t eps2= %f\n", eH, eO);
  PetscPrintf(PETSC_COMM_WORLD,"sig1= %f \t sig2= %f\n", sH, sO);
  PetscPrintf(PETSC_COMM_WORLD,"q1  = %f \t q2  = %f\n", qH, qO);
  PetscPrintf(PETSC_COMM_WORLD,"--------------------------------\n");

  BHD->rhos[0] = PD->rho; //2.*PD->rho;
  BHD->rhos[1] = PD->rho;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &BHD->fft_mat, &BHD->da, &BHD->dc);

  const DA da = BHD->da;         /* shorter alias */

  /* Create global vectors */
  DACreateGlobalVector(da, &BHD->g_ini[0]);
  DACreateGlobalVector(da, &BHD->g_ini[1]);
  DACreateGlobalVector(da, &BHD->gHO_ini);
  DACreateGlobalVector(da, &BHD->u2[0][0]);
  DACreateGlobalVector(da, &BHD->u2[1][1]);
  DACreateGlobalVector(da, &BHD->u2[0][1]);
  DACreateGlobalVector(da, &BHD->cH);
  DACreateGlobalVector(da, &BHD->cO);
  DACreateGlobalVector(da, &BHD->cHO);
  FOR_DIM
    {
      DACreateGlobalVector(da, &BHD->F[0][0][dim]);
      DACreateGlobalVector(da, &BHD->F[1][1][dim]);
      DACreateGlobalVector(da, &BHD->F[0][1][dim]);
      DACreateGlobalVector(da, &BHD->F_l[0][0][dim]);
      DACreateGlobalVector(da, &BHD->F_l[1][1][dim]);
      DACreateGlobalVector(da, &BHD->F_l[0][1][dim]);
      DACreateGlobalVector(da, &BHD->v[dim]);
    }

  /* Allocate memory for fft */
  FOR_DIM
    DACreateGlobalVector (BHD->dc, &BHD->fg2_fft[dim]);

  /* Complex scratch vector. FIXME: is it used in pure code? */
  DACreateGlobalVector (BHD->dc, &BHD->fft_scratch);
  DACreateGlobalVector (BHD->dc, &BHD->gfg2_fft);

  /* FIXME: these probably differ only by factors: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD->dc, &BHD->u2_fft[i][j]);
        BHD->u2_fft[j][i] = BHD->u2_fft[i][j];
      }

  return BHD;
}



static void BGY3dH2OData_free(State *BHD)
{
  MPI_Barrier( PETSC_COMM_WORLD);

  FOR_DIM
    {
      VecDestroy(BHD->F[0][0][dim]);
      VecDestroy(BHD->F[1][1][dim]);
      VecDestroy(BHD->F[0][1][dim]);
      VecDestroy(BHD->F_l[0][0][dim]);
      VecDestroy(BHD->F_l[1][1][dim]);
      VecDestroy(BHD->F_l[0][1][dim]);
      VecDestroy(BHD->v[dim]);
      VecDestroy (BHD->fg2_fft[dim]);
    }
  VecDestroy (BHD->gfg2_fft);

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      VecDestroy (BHD->u2_fft[i][j]);

  VecDestroy(BHD->g_ini[0]);
  VecDestroy(BHD->g_ini[1]);
  VecDestroy(BHD->gHO_ini);
  VecDestroy(BHD->u2[0][0]);
  VecDestroy(BHD->u2[1][1]);
  VecDestroy(BHD->u2[0][1]);
  VecDestroy(BHD->cH);
  VecDestroy(BHD->cO);
  VecDestroy(BHD->cHO);
  VecDestroy (BHD->fft_scratch);

#ifdef L_BOUNDARY
  MatDestroy (BHD->M);
  KSPDestroy(BHD->ksp);
#endif

  DADestroy (BHD->da);
  DADestroy (BHD->dc);
  MatDestroy (BHD->fft_mat);

  free(BHD);
}

static real LJ_repulsive(real r, real epsilon, real sigma)
{
  real sr6, sr, re;

  r=r+SHIFT;

  sr = sigma/r;
  sr6 = SQR(sr)*SQR(sr)*SQR(sr);

  re= 4.*epsilon*sr6*sr6;

  if(fabs(re)>epsilon*CUTOFF)
    return epsilon*CUTOFF;
  else
    return re;
}


/* NOTE: so far  in all cases the returned  result contains the factor
   q2. */
real Coulomb_short (real r, real q2)
{
    real re;

    if (r==0) {
        /* FIXME: pointless branch here: */
        if (q2 > 0)
            return EPSILON0INV * q2 * (CUTOFF*1.0e-5);
        else
            return EPSILON0INV * q2 * (CUTOFF*1.0e-5); //1.0e+4;
    }

    re = EPSILON0INV * q2 * erfc(G * r) / r;

    /* Check for large values remember: exp(-re) will be computed: */
    if (fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
        return EPSILON0INV * q2 * (CUTOFF*1.0e-5);
    else
        return re;
}

real Coulomb_short_grad( real r, real rx, real q2 )
{
  real re;

  if(rx==0)
    return 0;
  if(r==0)
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);


  re = - EPSILON0INV * q2 * (erfc(G*r) + 2.*G/sqrt(M_PI)*r*exp(-G*G*r*r))*rx/pow(r,3.0);

  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);
  else
    return re;
}

/* Coulomb_long  (r, 1.0) =  EPSILON0INV *  erf (G  * r)  / r  in most
   general  case.   An  extra  argument  q2  is  the  overall  scaling
   factor. */
real Coulomb_long (real r, real q2)
{
  real re;

   if (r == 0.0)
     {
       return EPSILON0INV * q2 * G * 2.0 / sqrt(M_PI);
     }

  re = EPSILON0INV * q2 * erf (G * r) / r;

  /* Check for large values, remember: exp(-re) will be computed */
  if (fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
    return EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
  else
    return re;
}

real Coulomb_long_grad( real r, real rx, real q2)
{
  real re;

  if(r==0)
    return 0;


  re = - EPSILON0INV * q2 * (erf(G*r)
                             - 2.*G/sqrt(M_PI)*r*exp(-G*G*r*r))*rx/pow(r,3.0);

  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);
  else
    return re;
}


real Coulomb( real r, real q2)
{
  real  re;

   if(r==0)
     {

       return EPSILON0INV * q2 * (CUTOFF*1.0e-5);


     }


  re = EPSILON0INV * q2 /r;

  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return EPSILON0INV * q2 * (CUTOFF*1.0e-5);
/*   else if( -re>1e+1) */
/*     return -1e+1; */
  else
    return re;
}

real Coulomb_grad( real r, real rx, real q2)
{
  real re;

  if( rx==0 )
    return 0;
  if(r==0)
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);

  re = - EPSILON0INV * q2 * rx/pow(r,3.0);



  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)) )
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);
  else
    return re;
}

/* g := exp[-(g0 + dg)], with a sanity check: */
void ComputeH2O_g (Vec g, Vec g0, Vec dg)
{
  int local_size;
  PetscScalar *g_vec, *dg_vec, *g0_vec;

  VecGetArray (g, &g_vec);
  VecGetArray (g0, &g0_vec);
  VecGetArray (dg, &dg_vec);

  VecGetLocalSize (g, &local_size);

  for (int i = 0; i < local_size; i++)
    {
      real k = - g0_vec[i] - dg_vec[i];
      g_vec[i] = exp(k);

      assert(!isinf(g_vec[i]) && !isnan(g_vec[i]));
    }

  VecRestoreArray (g, &g_vec);
  VecRestoreArray (g, &g0_vec);
  VecRestoreArray (dg, &dg_vec);
}

void Smooth_Function(State *BHD, Vec g, real RL, real RR, real shift)
{
  DA da;
  int x[3], n[3], i[3];
  PetscScalar ***g_vec;
  real r[3], r_s, h[3], interval[2], s, r_rl_2, rr_rl_2r, rr_rl_3;

  const ProblemData *PD = BHD->PD;
  da = BHD->da;

  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);


  DAVecGetArray(da, g, &g_vec);
   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {

/*        if(i[2]==0 || i[2]==n[2]-1 || i[1]==0 || i[1]==n[1]-1 ||  */
/*           i[0]==0 || i[0]==n[0]-1) */
/*          g_vec[i[2]][i[1]][i[0]]=0; */

          FOR_DIM
            r[dim] = i[dim]*h[dim]+interval[0];
          r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

          if( r_s > RR)
            s=0;
          else if ( r_s > RL)
            {
              r_rl_2 = SQR(r_s-RL);
              rr_rl_2r = (3.*RR-RL-2.*r_s);
              rr_rl_3 = pow(RR-RL,3);

              s= 1. - r_rl_2*rr_rl_2r/rr_rl_3;

            }
          else
            s=1;
          g_vec[i[2]][i[1]][i[0]] = s* (g_vec[i[2]][i[1]][i[0]]-shift) + shift;
        }
  DAVecRestoreArray(da, g, &g_vec);
}


/*
  Long range  pair potential Vec uc  is intent(out) here,  same as its
  FFT  transform uc_fft.  Vec  fc[] is  filled with  the corresponding
  force derived by means of FFT from the potential uc.

  No side effects.
*/
void ComputeFFTfromCoulomb (State *BHD,
                            Vec uc, Vec fc[3], /* intent(out) */
                            Vec uc_fft,    /* complex, intent(out) */
                            Vec fc_fft[3], /* complex, intent(out) */
                            real factor)
{
  int x[3], n[3], i[3], ic[3];

  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real L = PD->interval[1] - PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***uc_fft_, ***fc_fft_[3];
  DAVecGetArray (BHD->dc, uc_fft, &uc_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fc_fft[dim], &fc_fft_[dim]);

   /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          FOR_DIM
            {
              if (i[dim] <= N[dim] / 2)
                ic[dim] = i[dim];
              else
                ic[dim] = i[dim] - N[dim];
            }

          if (ic[0] == 0 && ic[1] == 0 && ic[2] == 0)
            {
              uc_fft_[i[2]][i[1]][i[0]] = 0.0; /* complex */
              FOR_DIM
                fc_fft_[dim][i[2]][i[1]][i[0]] = 0; /* complex */
            }
          else
            {
              const real k2 = (SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0])) / SQR(L);
              const real fac = EPSILON0INV / M_PI / k2;
              const real sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);

              /* Potential, complex with zero imaginary part: */
              uc_fft_[i[2]][i[1]][i[0]] = factor * sign * fac * exp(- k2 * SQR(M_PI) / SQR(G));

              /* Force, imaginary: */
              FOR_DIM
                fc_fft_[dim][i[2]][i[1]][i[0]] = 2.0 * M_PI * ic[dim] / L * (I * uc_fft_[i[2]][i[1]][i[0]]);
            }
        }
  DAVecRestoreArray (BHD->dc, uc_fft, &uc_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fc_fft[dim], &fc_fft_[dim]);

  /* FFT potential */
  MatMultTranspose (BHD->fft_mat, uc_fft, uc);
  VecScale (uc, 1./L/L/L);

  /* FFT force, seems to be always requested: */
  assert (fc != NULL);
  FOR_DIM
    {
      MatMultTranspose (BHD->fft_mat, fc_fft[dim], fc[dim]);
      VecScale (fc[dim], 1./L/L/L);
    }
}


static void RecomputeInitialData (State *BHD, real damp, real damp_LJ)
{
  DA da;
  PetscScalar ***gHini_vec, ***gOini_vec, ***gHOini_vec;
  PetscScalar ***(fH_vec[3]),***(fO_vec[3]),***(fHO_vec[3]);
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]),***(fHOl_vec[3]);
  PetscScalar ***cH_vec, ***cHO_vec, ***cO_vec;
  real r[3], r_s, h[3], interval[2], beta, L;
  int x[3], n[3], i[3], k;
  // local LJ params and charge product
  real q2H, q2O, q2HO;
  real epsilonH, epsilonO, epsilonHO;
  real sigmaH, sigmaO, sigmaHO;


  epsilonH = BHD->LJ_paramsH[0];
  epsilonO = BHD->LJ_paramsO[0];
  epsilonHO = BHD->LJ_paramsHO[0];

  sigmaH = BHD->LJ_paramsH[1];
  sigmaO = BHD->LJ_paramsO[1];
  sigmaHO = BHD->LJ_paramsHO[1];

  q2H = BHD->LJ_paramsH[2];
  q2O = BHD->LJ_paramsO[2];
  q2HO = BHD->LJ_paramsHO[2];


  const ProblemData *PD = BHD->PD;
  da = BHD->da;

  PetscPrintf(PETSC_COMM_WORLD,"Recomputing initial data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);


  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;

  real periodic[27][3]={{0,0,0},{L,0,0},{-L,0,0},{0,L,0},{0,-L,0},{0,0,L},{0,0,-L},
                       {L,L,0},{-L,L,0},{L,-L,0},{-L,-L,0},
                       {L,L,L},{-L,L,L},{L,-L,L},{-L,-L,L},
                       {L,L,-L},{-L,L,-L},{L,-L,-L},{-L,-L,-L},
                       {0,L,L},{0,-L,L},{0,L,-L},{0,-L,-L},
                       {L,0,L},{-L,0,L},{L,0,-L},{-L,0,-L}};



  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);


   FOR_DIM
    {
      VecSet(BHD->F[0][0][dim],0.0);
      VecSet(BHD->F[1][1][dim],0.0);
      VecSet(BHD->F[0][1][dim],0.0);
      VecSet(BHD->F_l[0][0][dim],0.0);
      VecSet(BHD->F_l[1][1][dim],0.0);
      VecSet(BHD->F_l[0][1][dim],0.0);
    }
  VecSet(BHD->g_ini[0], 0.0);
  VecSet(BHD->g_ini[1], 0.0);
  VecSet(BHD->gHO_ini, 0.0);

  /*********************************************/
  /* Compute fft from Coulomb potential (long) */
  /********************************************/

  /* FFT image of Coulomb long  force in BHD->fg2_fft[3] appears to be
     ignored since overwritten. So here these are work arrays: */
  ComputeFFTfromCoulomb (BHD, BHD->u2[0][1], BHD->F_l[0][1],
                         BHD->u2_fft[0][1], BHD->fg2_fft,
                         q2HO * damp);
  ComputeFFTfromCoulomb (BHD, BHD->u2[0][0], BHD->F_l[0][0],
                         BHD->u2_fft[0][0], BHD->fg2_fft,
                         q2H * damp);
  ComputeFFTfromCoulomb (BHD, BHD->u2[1][1], BHD->F_l[1][1],
                         BHD->u2_fft[1][1], BHD->fg2_fft,
                         q2O * damp);

  FOR_DIM
    {
      VecAXPY(BHD->F[0][1][dim], 1.0, BHD->F_l[0][1][dim]);
      VecAXPY(BHD->F[0][0][dim], 1.0, BHD->F_l[0][0][dim]);
      VecAXPY(BHD->F[1][1][dim], 1.0, BHD->F_l[1][1][dim]);
    }

  /**********************************************/

  DAVecGetArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecGetArray(da, BHD->g_ini[1], &gOini_vec);
  DAVecGetArray(da, BHD->gHO_ini, &gHOini_vec);
  DAVecGetArray(da, BHD->cH, &cH_vec);
  DAVecGetArray(da, BHD->cO, &cO_vec);
  DAVecGetArray(da, BHD->cHO, &cHO_vec);
  FOR_DIM
    {
      DAVecGetArray(da, BHD->F[0][0][dim], &fH_vec[dim]);
      DAVecGetArray(da, BHD->F[1][1][dim], &fO_vec[dim]);
      DAVecGetArray(da, BHD->F[0][1][dim], &fHO_vec[dim]);
      DAVecGetArray(da, BHD->F_l[0][0][dim], &fHl_vec[dim]);
      DAVecGetArray(da, BHD->F_l[1][1][dim], &fOl_vec[dim]);
      DAVecGetArray(da, BHD->F_l[0][1][dim], &fHOl_vec[dim]);
    }



  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          /* set force vectors */

          for(k=0; k<1; k++)
            {

          FOR_DIM
            r[dim] = i[dim]*h[dim]+interval[0]+periodic[k][dim];


          r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );


          /* Lennard-Jones */
          gHini_vec[i[2]][i[1]][i[0]] +=
            damp_LJ * beta* Lennard_Jones( r_s, epsilonH, sigmaH);
          gOini_vec[i[2]][i[1]][i[0]] +=
            damp_LJ * beta* Lennard_Jones( r_s, epsilonO, sigmaO);
          gHOini_vec[i[2]][i[1]][i[0]] +=
            damp_LJ * beta* Lennard_Jones( r_s, epsilonHO, sigmaHO);

          /* Coulomb short */
          gHini_vec[i[2]][i[1]][i[0]] +=
            damp*beta* Coulomb_short( r_s, q2H);
          gOini_vec[i[2]][i[1]][i[0]] +=
            damp*beta* Coulomb_short( r_s, q2O);
          gHOini_vec[i[2]][i[1]][i[0]] +=
            damp*beta* Coulomb_short( r_s, q2HO);

          /* Coulomb long */
/*        gHini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsH); */
/*        gOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */
/*        gHOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */

           /* Coulomb */
/*        gHini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsH); */
/*        gOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsO); */
/*        gHOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsHO); */

           FOR_DIM
            {
              /* Lennard-Jones */
              fH_vec[dim][i[2]][i[1]][i[0]] +=
                damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonH, sigmaH);
              fO_vec[dim][i[2]][i[1]][i[0]] +=
                damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonO, sigmaO);
              fHO_vec[dim][i[2]][i[1]][i[0]] +=
                damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonHO, sigmaHO);

              /* Coulomb short */
              fH_vec[dim][i[2]][i[1]][i[0]] +=
                damp * Coulomb_short_grad( r_s, r[dim], q2H);
              fO_vec[dim][i[2]][i[1]][i[0]] +=
                damp * Coulomb_short_grad( r_s, r[dim], q2O);
              fHO_vec[dim][i[2]][i[1]][i[0]] +=
                damp * Coulomb_short_grad( r_s, r[dim], q2HO);

              /* Coulomb long */
/*            fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/*            fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/*            fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

              /* Coulomb long */
/*            fHl_vec[dim][i[2]][i[1]][i[0]] =  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/*            fOl_vec[dim][i[2]][i[1]][i[0]] =  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/*            fHOl_vec[dim][i[2]][i[1]][i[0]] =  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

              /* Coulomb */
/*            fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsH); */
/*            fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsO); */
/*            fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsHO); */

              /* deterministic correction */
              cH_vec[i[2]][i[1]][i[0]] =
                exp(-beta* LJ_repulsive( r_s, epsilonH, sigmaH));
              cO_vec[i[2]][i[1]][i[0]] =
                exp(-beta* LJ_repulsive( r_s, epsilonO, sigmaO));
              cHO_vec[i[2]][i[1]][i[0]] =
                exp(-beta* LJ_repulsive( r_s, epsilonHO, sigmaHO));
            }
            }
        }

  DAVecRestoreArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecRestoreArray(da, BHD->g_ini[1], &gOini_vec);
  DAVecRestoreArray(da, BHD->gHO_ini, &gHOini_vec);
  DAVecRestoreArray(da, BHD->cH, &cH_vec);
  DAVecRestoreArray(da, BHD->cO, &cO_vec);
  DAVecRestoreArray(da, BHD->cHO, &cHO_vec);
  FOR_DIM
    {
      DAVecRestoreArray(da, BHD->F[0][0][dim], &fH_vec[dim]);
      DAVecRestoreArray(da, BHD->F[1][1][dim], &fO_vec[dim]);
      DAVecRestoreArray(da, BHD->F[0][1][dim], &fHO_vec[dim]);
      DAVecRestoreArray(da, BHD->F_l[0][0][dim], &fHl_vec[dim]);
      DAVecRestoreArray(da, BHD->F_l[1][1][dim], &fOl_vec[dim]);
      DAVecRestoreArray(da, BHD->F_l[0][1][dim], &fHOl_vec[dim]);
    }
}

/* All vectors are complex here. No side effects. */
static void kapply (const State *BHD,
                    Vec fg2_fft[3],          /* intent(in) */
                    Vec g_fft, Vec coul_fft, /* intent(in) */
                    Vec dg_fft)              /* intent(out) */
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real h = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];
  const real fac = L / (2. * M_PI); /* BHD->f ist nur grad U, nicht F=-grad U  */

  complex ***g_fft_, ***dg_fft_, ***coul_fft_, ***fg2_fft_[3];
  DAVecGetArray (BHD->dc, g_fft, &g_fft_);
  DAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecGetArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* Get local portion of the grid */
  int x[3], n[3], i[3], ic[3];
  DAGetCorners(BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          dg_fft_[i[2]][i[1]][i[0]] = 0.0; /* complex zero */

          FOR_DIM
            {
              if (i[dim] <= N[dim] / 2)
                ic[dim] = i[dim];
              else
                ic[dim] = i[dim] - N[dim];
            }

          if( ic[0]==0 && ic[1]==0 && ic[2]==0)
            {
              dg_fft_[i[2]][i[1]][i[0]] = 0.0; // coul1_fft[0].re*h*(g_fft[0].re-N[0]*N[1]*N[2]);
            }
          else
            {
              const real k = SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0]);
              const real k_fac = h * h * fac / k;

              /* phase shift factor for x=x+L/2 */
              const real sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);

              /* "I" is an imaginary unit here: */
              FOR_DIM
                dg_fft_[i[2]][i[1]][i[0]] += ic[dim] * k_fac * sign *
                (-I * fg2_fft_[dim][i[2]][i[1]][i[0]] * g_fft_[i[2]][i[1]][i[0]]);

              /* Long  range  Coulomb part.  Note  there  is not  "-I"
                 factor here: */
              dg_fft_[i[2]][i[1]][i[0]] += h * sign *
                (coul_fft_[i[2]][i[1]][i[0]] * g_fft_[i[2]][i[1]][i[0]]);
            }
        }
  DAVecRestoreArray (BHD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecRestoreArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);
}

/*
  Side effects:  uses BHD->{fft_scratch, fg2_fft[],  gfg2_fft} as work
  Vecs. Does 2 * (4 + 1) FFTs. Two inverse ones
*/
static void Compute_dg_H2O_inter (State *BHD,
                                  Vec f1[3], Vec f1_l[3], Vec g1a, Vec g1b,
                                  Vec coul1_fft, real rho1,
                                  Vec f2[3], Vec f2_l[3], Vec g2a, Vec g2b,
                                  Vec coul2_fft, real rho2,
                                  Vec dg, Vec dg_help)
{
  const ProblemData *PD = BHD->PD;
  const real L = PD->interval[1] - PD->interval[0];

  Vec g_fft = BHD->fft_scratch;
  Vec *fg2_fft = BHD->fg2_fft;  /* fg2_fft[3] */
  Vec dg_fft = BHD->gfg2_fft;

  /************************************************/
  /* rho1 * F1*g1a g1b */
  /************************************************/

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult (BHD->v[dim], g1a, f1[dim]);

      /* special treatment: Coulomb long */
      VecAXPY (BHD->v[dim], -1.0, f1_l[dim]);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }

  /* fft(g) */
  MatMult (BHD->fft_mat, g1b, g_fft);

  kapply (BHD, fg2_fft, g_fft, coul1_fft, dg_fft);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg_help);

  VecScale (dg_help, rho1 * PD->beta/L/L/L);

  VecCopy (dg_help, dg);

  /************************************************/
  /* rho2 * F2*g2a g2b */
  /************************************************/

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult (BHD->v[dim], g2a, f2[dim]);

      /* special treatment: Coulomb long */
      VecAXPY (BHD->v[dim], -1.0, f2_l[dim]);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }

  /* fft(g-1) */
  MatMult (BHD->fft_mat, g2b, g_fft);

  kapply (BHD, fg2_fft, g_fft, coul2_fft, dg_fft);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg_help);

  VecScale (dg_help, rho2 * PD->beta/L/L/L);

  VecAXPY (dg, 1.0, dg_help);
}


/* Compute  intramolecular  part. */
static void Compute_dg_H2O_intraIII (State *BHD, Vec f[3], Vec f_l[3], Vec g1, Vec tg,
                                     Vec coul_fft, real rab, Vec dg, Vec dg_help)
{
  int x[3], n[3], i[3], ic[3], local_size;
  real k_fac, k;
  PetscScalar *v_vec, *tg_vec;

  const ProblemData *PD = BHD->PD;
  const DA da = BHD->da;
  const int *N = PD->N;         /* N[3] */
  Vec *fg2_fft = BHD->fg2_fft;  /* fg2_fft[3] */

  const real h = PD->h[0] * PD->h[1] * PD->h[2];
  Vec dg_fft = BHD->gfg2_fft;
  const real L = PD->interval[1] - PD->interval[0];
  const real fac = L / (2. * M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/

  /* fft(f1*g1) */
  FOR_DIM
    {
      VecPointwiseMult(BHD->v[dim], g1, f[dim]);
      /* special treatment: Coulomb long */
      VecAXPY(BHD->v[dim], -1.0, f_l[dim]);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }


  struct {PetscScalar re, im;} ***fg2_fft_[3];
  FOR_DIM
    DAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* Compute int(F*g*omega) */
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {

          FOR_DIM
            {
              if( i[dim] <= N[dim]/2)
                ic[dim] = i[dim];
              else
                ic[dim] = i[dim] - N[dim];
            }

          if( ic[0]==0 && ic[1]==0 && ic[2]==0)
            {
              FOR_DIM
                fg2_fft_[dim][i[2]][i[1]][i[0]].re =0 ;
              FOR_DIM
                fg2_fft_[dim][i[2]][i[1]][i[0]].im =0 ;

            }
          else
            {
              k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
              k = 2.0*M_PI*sqrt(k)*rab/L;

              FOR_DIM
                fg2_fft_[dim][i[2]][i[1]][i[0]].re *= h * sin(k) / k;
              FOR_DIM
                fg2_fft_[dim][i[2]][i[1]][i[0]].im *= h * sin(k) / k;
            }
        }
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  FOR_DIM
    {
      MatMultTranspose (BHD->fft_mat, fg2_fft[dim], BHD->v[dim]);
      VecScale(BHD->v[dim], 1./L/L/L);
    }




  /* int(..)/tg */
  /* Laplace^-1 * divergence */


  /*******************************************/
  /*    int/tg  */
  /*******************************************/
  VecGetArray( tg, &tg_vec);
  FOR_DIM
    {
      VecGetLocalSize(BHD->v[dim], &local_size);
      VecGetArray( BHD->v[dim] , &v_vec);


      for (int index=0; index < local_size; index++)
        {
          k = tg_vec[index];
          if( k<NORM_REG)
            v_vec[index] = v_vec[index] / NORM_REG;
          else v_vec[index] = v_vec[index] / k;

/*        if( v_vec[index]<0) */
/*          { */
/*            if( fabs(k) <= 0.8) */
/*              v_vec[index] = v_vec[index] / 0.8; */
/*            else */
/*              v_vec[index] = v_vec[index] / k;         */

/*          } */

/*        else */
/*          { */
/*            if( fabs(k) <= 1.0e-4) */
/*              v_vec[index] = v_vec[index] / 1.0e-4 ; */
/*            else */
/*              v_vec[index] = v_vec[index] / k; */
/*            if( v_vec[index] > 20.) */
/*              v_vec[index]=20.0; */
/*          } */

        }
      VecRestoreArray( BHD->v[dim], &v_vec);

      //VecPointwiseDivide( BHD->v[dim], BHD->v[dim], tg);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }
  VecRestoreArray (tg, &tg_vec);

  /*******************************************/
  struct {PetscScalar re, im;} ***dg_fft_, ***coul_fft_;
  DAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecGetArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

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
              dg_fft_[i[2]][i[1]][i[0]].re = 0;
              dg_fft_[i[2]][i[1]][i[0]].im = 0;
              FOR_DIM
                {
                  fg2_fft_[dim][i[2]][i[1]][i[0]].im=0;
                  fg2_fft_[dim][i[2]][i[1]][i[0]].re=0;
                }
            }
          else
            {
              k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
              k_fac = fac/k;
              k = 2.0*M_PI*sqrt(k)*rab/L;

              FOR_DIM
                dg_fft_[i[2]][i[1]][i[0]].re += ic[dim] * k_fac * h * fg2_fft_[dim][i[2]][i[1]][i[0]].im;

              FOR_DIM
                dg_fft_[i[2]][i[1]][i[0]].im -= ic[dim] * k_fac * h * fg2_fft_[dim][i[2]][i[1]][i[0]].re;


              FOR_DIM
                fg2_fft_[dim][i[2]][i[1]][i[0]].re = -ic[dim] / fac * coul_fft_[i[2]][i[1]][i[0]].im * sin(k) / k;
              FOR_DIM
                fg2_fft_[dim][i[2]][i[1]][i[0]].im =  ic[dim] / fac * coul_fft_[i[2]][i[1]][i[0]].re * sin(k) / k;;
            }
        }
  DAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecRestoreArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg_help);

  VecScale(dg_help, PD->beta/L/L/L);

  /* back transformation of coulomb part, divide by tg and forward transfromation */
  VecGetArray( tg, &tg_vec);
  FOR_DIM
    {
      MatMultTranspose (BHD->fft_mat, fg2_fft[dim], BHD->v[dim]);
      VecScale(BHD->v[dim], 1./L/L/L);

     VecGetLocalSize(BHD->v[dim], &local_size);
     VecGetArray( BHD->v[dim] , &v_vec);


      for (int index=0; index < local_size; index++)
        {

          k = tg_vec[index];
          if( k<NORM_REG)
            v_vec[index] = v_vec[index] / NORM_REG;
          else v_vec[index] = v_vec[index] / k;
/*        if( v_vec[index]<0) */
/*          { */
/*            if( fabs(k) <= 0.8) */
/*              v_vec[index] = v_vec[index] / 0.8 ; */
/*            else */
/*              v_vec[index] = v_vec[index] / k; */
/*          } */
/*        else */
/*          { */
/*            if( fabs(k) <= 1.0e-4) */
/*              v_vec[index] = v_vec[index] / 1.0e-4 ; */
/*            else */
/*              v_vec[index] = v_vec[index] / k; */
/*            if( v_vec[index] > 20.) */
/*              v_vec[index]=20.0; */
/*          } */

        }
      VecRestoreArray( BHD->v[dim], &v_vec);

      //VecPointwiseDivide(BHD->v[dim], BHD->v[dim], tg);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }
  VecRestoreArray( tg, &tg_vec);

  DAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

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
              dg_fft_[i[2]][i[1]][i[0]].re = 0;
              dg_fft_[i[2]][i[1]][i[0]].im = 0;

            }
          else
            {
              k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
              k_fac = fac/k;

              FOR_DIM
                dg_fft_[i[2]][i[1]][i[0]].re += ic[dim] * k_fac * h * fg2_fft_[dim][i[2]][i[1]][i[0]].im;

              FOR_DIM
                dg_fft_[i[2]][i[1]][i[0]].im -= ic[dim] * k_fac * h * fg2_fft_[dim][i[2]][i[1]][i[0]].re;
            }
        }
  DAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg);

  VecScale(dg, PD->beta/L/L/L);

  VecAXPY(dg, 1.0, dg_help);
}


/*
  Compute   intramolecular  part.    The  first   thing  it   does  is
  transforming g(x) -> g(k). The real space representation g(x) is not
  otherwise used.

  Vec g is intent(in).
  Vec dg is intent(out).

  Side effects: uses BHD->fft_scratch as temp Vec.

  FIXME: compare the code to normalization_intra().
*/
void Compute_dg_H2O_intra_ln (State *BHD, Vec g, real rab, Vec dg)
{
  /* g(x) -> g(k): */
  MatMult (BHD->fft_mat, g, BHD->fft_scratch);

  /*
    FIXME:  argument aliasing!   BHD->fft_scratch is  the  input alias
    work  array.  Its  contents  is  destroyed on  return.  Vec dg  is
    intent(out) here:
  */
  normalization_intra (BHD, BHD->fft_scratch, rab, BHD->fft_scratch, dg);

  /* ln(g) */
  {
    PetscScalar *g_vec;
    int local_size;

    VecGetArray (dg, &g_vec);
    VecGetLocalSize (dg, &local_size);

    for (int i = 0; i < local_size; i++)
      {
        real g_i = g_vec[i];
        assert (g_i >= 1.0e-8);  /* See normalization_intra() */

        g_vec[i] = -log(g_i);
      }
    VecRestoreArray(dg, &g_vec);
  }

  /* Ensure normalization condition int(u)=0 */
  /* VecSum (dg, &k); */
  /* VecShift (dg, -k/N[0]/N[1]/N[2]); */
}


/*
  Compute normalization condition.  Output  is the convolution of g(x)
  with ω(x)  = δ(|x| - r)  / 4πr². This function  takes momentum space
  g(k)  as  input.   The  convolution  result is  manipulated  in  the
  real-space  rep, unfortunately.  Though this  manipulation  is plain
  screening of the too small (negative) values.

  Vec dg is  intent(out).

  Needs a complex  work vector. As a matter of fact,  the work Vec and
  the  input  Vec  g_fft  may  be  aliased. Then  the  input  will  be
  destroyed, of course.

  FIXME: compare the code to Compute_dg_H2O_intra_ln().
 */
static void normalization_intra (const State *BHD,
                                 Vec g_fft, /* g(k), intent(in) */
                                 real rab,
                                 Vec work, /* complex work Vec */
                                 Vec dg)   /* ĝ(x), intent(out) */
{
  int x[3], n[3], i[3], ic[3];

  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */

  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***g_fft_, ***work_;
  DAVecGetArray (BHD->dc, g_fft, &g_fft_);
  DAVecGetArray (BHD->dc, work, &work_);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          FOR_DIM
            {
              if (i[dim] <= N[dim] / 2)
                ic[dim] = i[dim];
              else
                ic[dim] = i[dim] - N[dim];
            }

          /* Get g(k): */
          complex gk = g_fft_[i[2]][i[1]][i[0]];

          /* Scale by ω(k): */
          if (ic[0] == 0 && ic[1] == 0 && ic[2] == 0)
            gk *= h3;         /* sinc(0) = 1 */
          else
            {
              const real k2 = SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0]);
              const real k = 2.0 * M_PI * sqrt(k2) * rab / L;

              /* + hier richtig !! */
              gk *= h3 * (sin(k) / k);
            }

          /* Set ĝ(k): */
          work_[i[2]][i[1]][i[0]] = gk;
        }
  DAVecRestoreArray (BHD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BHD->dc, work, &work_);

  /* Inverse FFT, ĝ(k) -> ĝ(x): */
  MatMultTranspose (BHD->fft_mat, work, dg);

  VecScale (dg, 1./L/L/L);

  /* g >= 0 ?? */
  {
    PetscScalar *g_vec;
    int local_size;

    VecGetArray (dg, &g_vec);
    VecGetLocalSize (dg, &local_size);

    for (int i = 0; i < local_size; i++)
      {
        real g_i = g_vec[i];
        if (g_i < 1.0e-8)
          g_i = 1.0e-8;

        g_vec[i] = g_i;
      }
    VecRestoreArray (dg, &g_vec);
  }
}

/*
 * FIXME: consider using normalization_intra() instead.
 *
 * Compute  normalization  condition. This  function  appears to  take
 * g(x), FFT it into g(k) and operate with the latter exclusively. The
 * result is manipulated in  the real-space rep, unfortunately. Though
 * this manipulation  is plain screening  of the too  small (negative)
 * values.
 *
 * Vec dg is  intent(out).
 *
 * Side effects: uses BHD->fft_scratch as work Vec.
 */
static void Compute_dg_H2O_normalization_intra (const State *BHD, Vec g, real rab,
                                                Vec dg) /* intent(out) */
{
  /* fft(g/t) */
  MatMult (BHD->fft_mat, g, BHD->fft_scratch);

  /*
    FIXME: argument aliasing! BHD->fft_scratch is used as input and as
    work array here.  Gets  overwritten as the result.  Huge potential
    for confusion:
  */
  normalization_intra (BHD, BHD->fft_scratch, rab, BHD->fft_scratch,
                       dg); /* dg is intent(out) */
}


/* w := x / max(y, thresh), essentially: */
static void safe_pointwise_divide (Vec w, /* intent(out) */
                                   Vec x, /* intent(in) */
                                   Vec y, /* intent(in) */
                                   real thresh)
{
  int local_size;
  PetscScalar *w_vec, *x_vec, *y_vec;

  VecGetLocalSize (w, &local_size);

  VecGetArray (w, &w_vec);
  VecGetArray (x, &x_vec);
  VecGetArray (y, &y_vec);

  for (int i = 0; i < local_size; i++) {
      real y_i = y_vec[i];
      if (y_i < thresh)
          w_vec[i] = x_vec[i] / thresh;
      else
          w_vec[i] = x_vec[i] / y_i;
  }

  VecRestoreArray (w, &w_vec);
  VecRestoreArray (x, &x_vec);
  VecRestoreArray (y, &y_vec);
}


/*
  Vec t is intent(out).
  Vec dg is intent(out).
*/
static void Solve_NormalizationH2O_small (const State *BHD, Vec gc, real rc, Vec g,
                                          Vec t,  /* intent(out) */
                                          Vec dg) /* intent(out) */
{
  /* Vec dg, isintent(out) here: */
  Compute_dg_H2O_normalization_intra (BHD, gc, rc, dg);

  /* t(x)  = g(x)  / ĝ(x),  avoiding small  denominators. Some  of the
     commented code used VecPointwiseDivide(t, g, dg). */
  safe_pointwise_divide (t, g, dg, NORM_REG2);

  bgy3d_boundary_set (BHD, dg, 1.0);
}


void bgy3d_solve_normalization (const State *BHD,
                                Vec gc_fft, /* intent(in) */
                                real rc,
                                Vec g,     /* intent(in) */
                                Vec t)     /* intent(out) */
{
  /*
    ĝ(x) goes into Vec t  which is is intent(out) here, fft_scratch is
    a work array:
  */
  normalization_intra (BHD, gc_fft, rc, BHD->fft_scratch, t);

  /*
    t(x) =  g(x) / ĝ(x)  (or rather t(x)  = g(x) / t(x)  with argument
    aliasing) avoiding small denominators.  Some of the commented code
    used VecPointwiseDivide(t, g, dg).
  */
  safe_pointwise_divide (t, g, t, NORM_REG2); /* argument aliasing */
}


/* solve */
Vec BGY3d_solve_2site (const ProblemData *PD, Vec g_ini)
{
  State *BHD;
  Vec dg_new, dg_new2, f;
  Vec g0[2][2], g[2][2], dg[2][2]; /* pair distributions, ij = ji */
  Vec tH, tO, tHO, tOH;
  real a=0.9, damp, damp_LJ;
  real aH, aHO, aO, a1=0.5, count=0.0;
  int iter, mycount=0, upwards=0, namecount=0;
  char nameO[20], nameH[20], nameHO[20];
  PetscScalar dgH_norm, dgO_norm, dgHO_norm;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Fourier ansatz...\n");

  BHD = BGY3dH2OData_Pair_malloc (PD);

  if (r_HH > 0)
    PetscPrintf (PETSC_COMM_WORLD, "WARNING: Solvent not a 2-Site model!\n");

  /* Compute initial data */
  RecomputeInitialData (BHD, 1.0, 1.0);

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Mixing parameter: */
  const real a0 = PD->lambda;

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD->da, &g[i][j]);
        g[j][i] = g[i][j];

        DACreateGlobalVector(BHD->da, &dg[i][j]);
        dg[j][i] = dg[i][j];
      }

  DACreateGlobalVector(BHD->da, &dg_new);
  DACreateGlobalVector(BHD->da, &dg_new2);
  DACreateGlobalVector(BHD->da, &f);

  DACreateGlobalVector(BHD->da, &tH);
  DACreateGlobalVector(BHD->da, &tO);
  DACreateGlobalVector(BHD->da, &tHO);
  DACreateGlobalVector(BHD->da, &tOH);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD->da, BHD->PD, &BHD->M, &BHD->ksp);

  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem  across iterations. Though  by now  I am  not sure  if KSP
    solver  really uses  them as  initial guess  as opposed  to always
    starting iterations from zero:
  */
  Vec x_lapl[2][2];             /* real */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD->da, &x_lapl[i][j]);
        VecSet (x_lapl[i][j], 0.0);
        x_lapl[j][i] = x_lapl[i][j]; /* pairwise property */
      }
#endif

  g0[0][0] = BHD->g_ini[0];
  g0[1][1] = BHD->g_ini[1];
  g0[0][1] = BHD->gHO_ini;

  /* set initial guess*/
  VecSet(dg[0][0],0);
  VecSet(dg[1][1],0);
  VecSet(dg[0][1],0);

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O"))
    bgy3d_read_g2 (2, dg, "dg%d%d.bin");

  VecSet(dg_new,0.0);

  for( damp=damp_start; damp <=1; damp+=0.1)
    {
      if(damp==-0.01)
        {
          damp_LJ=0;
          RecomputeInitialData (BHD, 0, 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }
      else if(damp==0.0)
        {
          damp_LJ=1.0;
          RecomputeInitialData (BHD, 0, 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }
      else
        {
          damp_LJ=1.0;
          count+=1.0;
          RecomputeInitialData (BHD, (damp), 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }

      bgy3d_impose_laplace_boundary (BHD, g0[0][0], tH, x_lapl[0][0]);
      bgy3d_impose_laplace_boundary (BHD, g0[1][1], tH, x_lapl[1][1]);
      bgy3d_impose_laplace_boundary (BHD, g0[0][1], tH, x_lapl[0][1]);

      /* g=g0*exp(-dg) */
      ComputeH2O_g (g[0][1], g0[0][1], dg[0][1]);
      ComputeH2O_g (g[0][0], g0[0][0], dg[0][0]);
      ComputeH2O_g (g[1][1], g0[1][1], dg[1][1]);

      /* Not sure if 0.0 as inital value is right. */
      real dgH_old = 0.0;
      real dgO_old = 0.0;
      real dgHO_old = 0.0;

      a1=a0;
      a=a0;
      aH=a;
      aO=a;
      aHO=a;

  for(iter=0; iter<max_iter; iter++)
    {



      if( !(iter%10) && iter>0 )
        {
          a=a1;
          aH=a;
          aO=a;
          aHO=a;
        }
      else if( iter==20)
        {
          aH=a;
          aO=a;
          aHO=a;
        }
      else
        {
          a=a0;
          aH=a;
          aO=a;
          aHO=a;
        }

      PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms: %e %e ", iter+1, NORM_REG, NORM_REG2);
      /* f=integral(g) */
      if (1)                    /* kflg was set with -pair */
        {
          /* g_OH */
          Compute_dg_H2O_inter(BHD,
                               BHD->F[1][1], BHD->F_l[1][1], g[1][1],
                               g[0][1],
                               BHD->u2_fft[1][1], BHD->rhos[1],
                               BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                               g[0][0],
                               BHD->u2_fft[0][1], BHD->rhos[0],
                               dg_new, f);
          VecScale(dg_new,damp_LJ);
          VecPointwiseMult(dg_new, dg_new, BHD->cHO);


          Solve_NormalizationH2O_small (BHD, g[0][1], r_HO, g[0][0], tH , dg_new2);
          Compute_dg_H2O_intra_ln(BHD, tH, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_small (BHD, g[0][1], r_HO, g[1][1], tO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[1][1], r_HO, tHO);
          Compute_dg_H2O_intraIII(BHD, BHD->F[1][1], BHD->F_l[1][1], tO, tHO,
                                 BHD->u2_fft[1][1], r_HO, dg_new2, f);
          VecAXPY(dg_new, 1.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[0][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, tH, x_lapl[0][1]);

          VecCopy(dg[0][1], f);
          VecAXPBY(dg[0][1], aHO, (1-aHO), dg_new);
          VecAXPY(f, -1.0, dg[0][1]);
          VecNorm(f, NORM_INFINITY, &dgHO_norm);
          PetscPrintf(PETSC_COMM_WORLD,"HO= %e  (%f)  ",  dgHO_norm/aHO, aHO);

          /* g_H */
          Compute_dg_H2O_inter(BHD,
                               BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                               g[0][1],
                               BHD->u2_fft[0][1], BHD->rhos[1],
                               BHD->F[0][0], BHD->F_l[0][0], g[0][0],
                               g[0][0],
                               BHD->u2_fft[0][0], BHD->rhos[0],
                               dg_new, f);
          VecScale(dg_new,damp_LJ);
          VecPointwiseMult(dg_new, dg_new, BHD->cH);

          Solve_NormalizationH2O_small (BHD, g[0][0], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_small (BHD, g[0][0], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[0][1], r_HO, tH);
          Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, tH,
                                 BHD->u2_fft[0][1], r_HO, dg_new2, f);
          VecAXPY(dg_new, 1.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[0][0]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, tH, x_lapl[0][0]);

          VecCopy(dg[0][0], f);
          VecAXPBY(dg[0][0], aH, (1-aH), dg_new);
          VecAXPY(f, -1.0, dg[0][0]);
          VecNorm(f, NORM_INFINITY, &dgH_norm);
          PetscPrintf(PETSC_COMM_WORLD,"H= %e  (%f)  ", dgH_norm/aH, aH);

          /* g_O */
          Compute_dg_H2O_inter(BHD,
                               BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                               g[0][1],
                               BHD->u2_fft[0][1], BHD->rhos[0],
                               BHD->F[1][1], BHD->F_l[1][1], g[1][1],
                               g[1][1],
                               BHD->u2_fft[1][1], BHD->rhos[1],
                               dg_new, f);
          VecScale(dg_new,damp_LJ);
          VecPointwiseMult(dg_new, dg_new, BHD->cO);

          Solve_NormalizationH2O_small (BHD, g[1][1], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_small (BHD, g[1][1], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[0][1], r_HO, tO);
          Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, tO,
                                 BHD->u2_fft[0][1], r_HO, dg_new2, f);
          VecAXPY(dg_new, 1.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[1][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, tH, x_lapl[1][1]);

          VecCopy(dg[1][1], f);
          VecAXPBY(dg[1][1], aO, (1-aO), dg_new);
          VecAXPY(f, -1.0,  dg[1][1]);
          VecNorm(f, NORM_INFINITY, &dgO_norm);
          PetscPrintf(PETSC_COMM_WORLD,"O= %e  (%f)  ", dgO_norm/aO, aO);

          /* ende: */
          ComputeH2O_g (g[0][1], g0[0][1], dg[0][1]);
          ComputeH2O_g (g[0][0], g0[0][0], dg[0][0]);
          ComputeH2O_g (g[1][1], g0[1][1], dg[1][1]);
        } /* of if (1) */

      /* (fancy) step size control */
      mycount++;
      if( ((iter-1)%10) &&
          (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
        {
          upwards=1;
        }
      else if(iter>20 && !((iter-1)%10) && upwards==0 &&
              (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
        {
          a1 /=2.;
          if(a1<a0)
            a1=a0;
          mycount=0;
        }
      else
        upwards=0;

      if(mycount>20 && a1<=0.5)
        {
          a1*=2.;
          mycount=0;
        }
      else if(mycount>20 && a1>0.5)
        {
          a1=1.0;
          mycount=0;
        }
      PetscPrintf(PETSC_COMM_WORLD,"count= %d  upwards= %d", mycount, upwards);
      dgH_old = dgH_norm/aH;
      dgO_old = dgO_norm/aO;
      dgHO_old = dgHO_norm/aHO;
      /*********************************/


      /* Last print didnt CR LF: */
      PetscPrintf (PETSC_COMM_WORLD, "\n");

      if (dgH_norm / aH<=norm_tol && dgO_norm / aO <= norm_tol && dgHO_norm / aHO <= norm_tol)
          break;
    }

  /* output */
  namecount++;
  sprintf(nameH, "vec00-%d.m", namecount-1);
  sprintf(nameO, "vec11-%d.m", namecount-1);
  sprintf(nameHO, "vec01-%d.m", namecount-1);

  PetscPrintf (PETSC_COMM_WORLD, "Writing files...");
  bgy3d_save_vec_ascii (nameH, g[0][0]); /* g_H */
  bgy3d_save_vec_ascii (nameO, g[1][1]); /* g_O */
  bgy3d_save_vec_ascii (nameHO, g[0][1]); /* g_HO */
  PetscPrintf(PETSC_COMM_WORLD,"done.\n");

  /* Save dg[][] to binary files: */
  if (bgy3d_getopt_test ("--save-H2O"))
    bgy3d_save_g2 (2, dg, "dg%d%d.bin");

  /* Save g2[][] to binary files: */
  bgy3d_save_g2 (2, g, "g%d%d.bin");

    }

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        VecDestroy (g[i][j]);
        VecDestroy (dg[i][j]);
        VecDestroy (x_lapl[i][j]);
      }

  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(f);

  VecDestroy(tH);
  VecDestroy(tO);
  VecDestroy(tHO);
  VecDestroy(tOH);

  BGY3dH2OData_free(BHD);

  return PETSC_NULL;
}


/* solve with product ansatz g=g0*dg */
Vec BGY3d_solve_3site (const ProblemData *PD, Vec g_ini)
{
  State *BHD;
  Vec dg_new, dg_new2, f;
  Vec g0[2][2], g[2][2], dg[2][2]; /* pair distributions, ij = ji */
  Vec tH, tO, tHO, tOH;
  real a=0.9, damp, damp_LJ;
  real aH, aHO, aO, a1=0.5, count=0.0;
  int iter, mycount=0, upwards=0, namecount=0;
  char nameO[20], nameH[20], nameHO[20];
  PetscScalar dgH_norm, dgO_norm, dgHO_norm;
  PetscScalar dgH_old, dgHO_old, dgO_old;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Fourier ansatz...\n");

  BHD = BGY3dH2OData_Pair_malloc (PD);

  if (r_HH < 0)
    {
      PetscPrintf (PETSC_COMM_WORLD, "Solvent not a 3-Site model!\n");
      exit(1);
    }

  /* Compute initial data */
  RecomputeInitialData (BHD, 1.0, 1.0);

  BHD->rhos[0] = 2.*BHD->rhos[0];

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Mixing parameter: */
  const real a0 = PD->lambda;

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD->da, &g[i][j]);
        g[j][i] = g[i][j];

        DACreateGlobalVector (BHD->da, &dg[i][j]);
        dg[j][i] = dg[i][j];
      }

  DACreateGlobalVector(BHD->da, &dg_new);
  DACreateGlobalVector(BHD->da, &dg_new2);
  DACreateGlobalVector(BHD->da, &f);

  DACreateGlobalVector(BHD->da, &tH);
  DACreateGlobalVector(BHD->da, &tO);
  DACreateGlobalVector(BHD->da, &tHO);
  DACreateGlobalVector(BHD->da, &tOH);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD->da, BHD->PD, &BHD->M, &BHD->ksp);

  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem  across iterations. Though  by now  I am  not sure  if KSP
    solver  really uses  them as  initial guess  as opposed  to always
    starting iterations from zero:
  */
  Vec x_lapl[2][2];             /* real */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD->da, &x_lapl[i][j]);
        VecSet (x_lapl[i][j], 0.0);
        x_lapl[j][i] = x_lapl[i][j]; /* pairwise property */
      }
#endif

  g0[0][0] = BHD->g_ini[0];
  g0[1][1] = BHD->g_ini[1];
  g0[0][1] = BHD->gHO_ini;

  /* set initial guess*/
  VecSet(dg[0][0],0);
  VecSet(dg[1][1],0);
  VecSet(dg[0][1],0);

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O"))
    bgy3d_read_g2 (2, dg, "dg%d%d.bin");

  VecSet(dg_new,0.0);

  for( damp=damp_start; damp <=damp_start; damp+=0.1)
    {
      if(damp==-0.01)
        {
          damp_LJ=0;
          RecomputeInitialData (BHD, 0, 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }
      else if(damp==0.0)
        {
          damp_LJ=1.0;
          RecomputeInitialData (BHD, 0, 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }
      else
        {
          damp_LJ=1.0;
          count+=1.0;
          RecomputeInitialData (BHD, (damp), 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }

      bgy3d_impose_laplace_boundary (BHD, g0[0][0], tH, x_lapl[0][0]);
      bgy3d_impose_laplace_boundary (BHD, g0[1][1], tH, x_lapl[1][1]);
      bgy3d_impose_laplace_boundary (BHD, g0[0][1], tH, x_lapl[0][1]);

      /* g=g0*exp(-dg) */
      ComputeH2O_g (g[0][1], g0[0][1], dg[0][1]);
      ComputeH2O_g (g[0][0], g0[0][0], dg[0][0]);
      ComputeH2O_g (g[1][1], g0[1][1], dg[1][1]);

      a1=a0;
      a=a0;
      aH=a;
      aO=a;
      aHO=a;

  for(iter=0; iter<max_iter; iter++)
    {
      if( !(iter%10) && iter>0 )
        {
          a=a1;
          aH=a;
          aO=a;
          aHO=a;
        }
      else if( iter==20)
        {
          aH=a;
          aO=a;
          aHO=a;
        }
      else
        {
          a=a0;
          aH=a;
          aO=a;
          aHO=a;
        }

      PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms:\t", iter+1);
      /* f=integral(g) */
      if (1)                    /* kflg was set when with -pair */
        {
          /* g_OH */
          Compute_dg_H2O_inter(BHD,
                               BHD->F[1][1], BHD->F_l[1][1], g[1][1],
                               g[0][1],
                               BHD->u2_fft[1][1], BHD->rhos[1],
                               BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                               g[0][0],
                               BHD->u2_fft[0][1], BHD->rhos[0],
                               dg_new, f);
          VecScale(dg_new,damp_LJ);

          Solve_NormalizationH2O_small (BHD, g[0][1], r_HO, g[0][0], tH , dg_new2);
          Compute_dg_H2O_intra_ln(BHD, tH, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 2.0, dg_new2);

          /* tO = gHO/int(gHO wHH) */
          Solve_NormalizationH2O_small (BHD, g[0][1], r_HH, g[0][1], tO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[0][1], r_HH, tHO);
          Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tO, tHO,
                                 BHD->u2_fft[0][1], r_HH, dg_new2, f);
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_small (BHD, g[0][1], r_HO, g[1][1], tO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[1][1], r_HO, tHO);
          Compute_dg_H2O_intraIII(BHD, BHD->F[1][1], BHD->F_l[1][1], tO, tHO,
                                 BHD->u2_fft[1][1], r_HO, dg_new2, f);
          VecAXPY(dg_new, 1.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[0][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, tH, x_lapl[0][1]);

          VecCopy(dg[0][1], f);
          VecAXPBY(dg[0][1], aHO, (1-aHO), dg_new);
          VecAXPY(f, -1.0, dg[0][1]);
          VecNorm(f, NORM_INFINITY, &dgHO_norm);
          PetscPrintf(PETSC_COMM_WORLD,"HO= %e  (%f)  ",  dgHO_norm/aHO, aHO);

          /* g_H */
          Compute_dg_H2O_inter(BHD,
                               BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                               g[0][1],
                               BHD->u2_fft[0][1], BHD->rhos[1],
                               BHD->F[0][0], BHD->F_l[0][0], g[0][0],
                               g[0][0],
                               BHD->u2_fft[0][0], BHD->rhos[0],
                               dg_new, f);
          VecScale(dg_new,damp_LJ);

          Solve_NormalizationH2O_small (BHD, g[0][0], r_HH, g[0][0], tH , dg_new2);
          Compute_dg_H2O_intra_ln(BHD, tH, r_HH, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_small (BHD, g[0][0], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

          /* tO = gH/int(gH wHH) */
          Solve_NormalizationH2O_small (BHD, g[0][0], r_HH, g[0][0], tO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[0][0], r_HH, tH);
          Compute_dg_H2O_intraIII(BHD, BHD->F[0][0], BHD->F_l[0][0], tO, tH,
                                  BHD->u2_fft[0][0], r_HH, dg_new2, f);
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_small (BHD, g[0][0], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[0][1], r_HO, tH);
          Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, tH,
                                 BHD->u2_fft[0][1], r_HO, dg_new2, f);
          VecAXPY(dg_new, 1.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[0][0]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, tH, x_lapl[0][0]);

          VecCopy(dg[0][0], f);
          VecAXPBY(dg[0][0], aH, (1-aH), dg_new);
          VecAXPY(f, -1.0, dg[0][0]);
          VecNorm(f, NORM_INFINITY, &dgH_norm);
          PetscPrintf(PETSC_COMM_WORLD,"H= %e  (%f)  ", dgH_norm/aH, aH);

          /* g_O */
          Compute_dg_H2O_inter(BHD,
                               BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                               g[0][1],
                               BHD->u2_fft[0][1], BHD->rhos[0],
                               BHD->F[1][1], BHD->F_l[1][1], g[1][1],
                               g[1][1],
                               BHD->u2_fft[1][1], BHD->rhos[1],
                               dg_new, f);
          VecScale(dg_new,damp_LJ);

          Solve_NormalizationH2O_small (BHD, g[1][1], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 2.0, dg_new2);

          Solve_NormalizationH2O_small (BHD, g[1][1], r_HO, g[0][1], tHO , dg_new2);
          Compute_dg_H2O_normalization_intra( BHD, g[0][1], r_HO, tO);
          Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, tO,
                                 BHD->u2_fft[0][1], r_HO, dg_new2, f);
          VecAXPY(dg_new, 2.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[1][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, tH, x_lapl[1][1]);

          VecCopy(dg[1][1], f);
          VecAXPBY(dg[1][1], aO, (1-aO), dg_new);
          VecAXPY(f, -1.0,  dg[1][1]);
          VecNorm(f, NORM_INFINITY, &dgO_norm);
          PetscPrintf(PETSC_COMM_WORLD,"O= %e  (%f)  ", dgO_norm/aO, aO);

          /* ende: */
          ComputeH2O_g (g[0][1], g0[0][1], dg[0][1]);
          ComputeH2O_g (g[0][0], g0[0][0], dg[0][0]);
          ComputeH2O_g (g[1][1], g0[1][1], dg[1][1]);
        } /* of if (1) */

      /* (fancy) step size control */
      mycount++;
      if( ((iter-1)%10) &&
          (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
        {
          upwards=1;
        }
      else if(iter>20 && !((iter-1)%10) && upwards==0 &&
              (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
        {
          a1 /=2.;
          if(a1<a0)
            a1=a0;
          mycount=0;
        }
      else
        upwards=0;

      if(mycount>20 && a1<=0.5)
        {
          a1*=2.;
          mycount=0;
        }
      else if(mycount>20 && a1>0.5)
        {
          a1=1.0;
          mycount=0;
        }
      PetscPrintf(PETSC_COMM_WORLD,"count= %d  upwards= %d", mycount, upwards);
      dgH_old = dgH_norm/aH;
      dgO_old = dgO_norm/aO;
      dgHO_old = dgHO_norm/aHO;
      /*********************************/

      if(dgH_norm/aH<=norm_tol &&  dgO_norm/aO<=norm_tol && dgHO_norm/aHO<=norm_tol)
        break;

      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }

  /* output */
  namecount++;
  sprintf(nameH, "vec00-%d.m", namecount-1);
  sprintf(nameO, "vec11-%d.m", namecount-1);
  sprintf(nameHO, "vec01-%d.m", namecount-1);

  PetscPrintf (PETSC_COMM_WORLD, "Writing files...");
  bgy3d_save_vec_ascii (nameH, g[0][0]);
  bgy3d_save_vec_ascii (nameO, g[1][1]);
  bgy3d_save_vec_ascii (nameHO, g[0][1]);
  PetscPrintf(PETSC_COMM_WORLD,"done.\n");

  /* Save dg[][] to binary files: */
  if (bgy3d_getopt_test ("--save-H2O"))
    bgy3d_save_g2 (2, dg, "dg%d%d.bin");

  /* Save g2[][] to binary files: */
  bgy3d_save_g2 (2, g, "g%d%d.bin");

    }

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        VecDestroy (g[i][j]);
        VecDestroy (dg[i][j]);
        VecDestroy (x_lapl[i][j]);
      }

  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(f);

  VecDestroy(tH);
  VecDestroy(tO);
  VecDestroy(tHO);
  VecDestroy(tOH);

  BGY3dH2OData_free(BHD);

  return PETSC_NULL;
}

