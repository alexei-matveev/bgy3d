/*==========================================================*/
/*  $Id: bgy3dH2O_solutes.c,v 1.3 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d_SolventParameters.h"

#define MAXATOM 20

typedef struct Solute
{
  char names[MAXATOM][5];       /* atom types */
  real x[MAXATOM][3];           /* the coordinates  */
  real sigma[MAXATOM];          /* sigma for LJ */
  real epsilon[MAXATOM];        /* epsilon for LJ */
  real q[MAXATOM];              /* charges */
  int max_atoms;
} Solute;

#if 0
static void ComputeSoluteDatafromCoulomb (BGY3dH2OData BHD, Vec uc, const real x0[3], real q2, real damp);
#endif
static void ComputeSoluteDatafromCoulombII (BGY3dH2OData BHD, Vec uc, const real x0[3], real q2, real damp);
static void ComputeSoluteDatafromCoulomb_QM (BGY3dH2OData BHD, Vec uc, Vec gs, real q, real damp);
static void RecomputeInitialSoluteData_QM (BGY3dH2OData BHD, const Solute *S, real damp, real damp_LJ);
static void RecomputeInitialSoluteData_II (BGY3dH2OData BHD, const Solute *S, real damp, real damp_LJ);

/*
 * These two functions  obey the same interface. They  are supposed to
 * get (1)  parameters of  the solvent site  such as its  location and
 * force field  parameters, and  (2) a description  of the  solute and
 * return a  real number such as  an interaction energy  or the charge
 * density:
 */
static real ljc (real x, real y, real z, real epsilon, real sigma, real charge, const Solute *S);
static real rho (real x, real y, real z, real epsilon, real sigma, real charge, const Solute *S);

/*
 * This function expects a callback obeying the above interface as one
 * of the arguments:
 */
static void field (DA da, const ProblemData *PD,
                   const Solute *S,
                   real epsilon, real sigma, real charge, real fact,
                   double (*f)(real x, real y, real z,
                               real eps, real sig, real chg,
                               const Solute *S),
                   Vec v);

/*********************************/
/* Water */
/*********************************/
static Solute Water =
  {
    { "O","OH","OH"},
    { {-0.2929, 0.000, 0.000},
      { 0.2929, 0.757, 0.000},
      { 0.2929,-0.757, 0.000}},
    {3.1506, 0.400, 0.400},
    {0.1521, 0.046, 0.046},
    {-0.834, 0.417, 0.417},
    3
  };


void RecomputeInitialSoluteData_Water(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Water.\n");
  RecomputeInitialSoluteData_II(BHD, &Water, damp, damp_LJ);
}

/*********************************/
/* CS2 */
/*********************************/

static Solute CarbonDisulfide =
  {
    { "C","S1","S2"},
    { {0.0, 0.0, 0.0},
      {-1.56, 0.0, 0.0},
      { 1.56, 0.0, 0.0}},
    {3.200, 3.520, 3.520},
    {0.10128, 0.39500, 0.39500},
    {-0.308, 0.154, 0.154},
    3
  };

void RecomputeInitialSoluteData_CS2(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is CarbonDisulfide.\n");
  RecomputeInitialSoluteData_II(BHD, &CarbonDisulfide, damp, damp_LJ);
}

/*********************************/
/* HCl */
/*********************************/

static Solute HydrogenChloride =
  {
    { "H","Cl"},
    { { 0.6285, 0.0, 0.0},
      {-0.6285, 0.0, 0.0}},
    {2.735, 3.353},
    {0.03971, 0.51434},
    {0.2, -0.2},
    2
  };

void RecomputeInitialSoluteData_HCl(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is HCl.\n");
#ifndef QM
  RecomputeInitialSoluteData_II(BHD, &HydrogenChloride, damp, damp_LJ);
#else
  RecomputeInitialSoluteData_QM(BHD, &HydrogenChloride, damp, damp_LJ);
#endif
}



/*********************************/
/* Methanol */
/*********************************/
/* static Solute Methanol =  */
/*   { */
/*     { "C","HC1","HC2","HC3","O","OH"},   */
/*     { {-0.748,-0.015, 0.024},  */
/*       {-1.293,-0.202,-0.901}, */
/*       {-1.263, 0.754, 0.600}, */
/*       {-0.699,-0.934, 0.609}, */
/*       { 0.558, 0.420,-0.278}, */
/*       { 0.716, 1.404, 0.137} }, */
/*     {3.6705, 2.3876, 2.3876, 2.3876, 3.1538, 0.4000},        */
/*     {0.0800, 0.0240, 0.0240, 0.0240, 0.1521, 0.0460},        */
/*     {-0.040, 0.090 , 0.090 , 0.090 , -0.660, 0.430 },         */
/*     6 */
/*   }; */

static Solute Methanol =
  {
    { "C","HC1","HC2","HC3","O","OH"},
    { {-0.748,-0.015, 0.024},
      {-1.293,-0.202,-0.901},
      {-1.263, 0.754, 0.600},
      {-0.699,-0.934, 0.609},
      { 0.558, 0.420,-0.278},
      { 0.716, 1.404, 0.137} },
/*     {3.500, 2.5000, 2.5000, 2.5000, 3.1200, 2.500},        */
/*     {0.066, 0.0300, 0.0300, 0.0300, 0.1700, 0.030},        */
/*     {0.145, 0.0400, 0.0400, 0.0400, -0.683, 0.418 },         */
    {3.500, 2.5000, 2.5000, 2.5000, 3.1200, 0.400},
    {0.066, 0.0300, 0.0300, 0.0300, 0.1700, 0.040},
    {0.145, 0.0400, 0.0400, 0.0400, -0.683, 0.418 },
    6
  };

/* static Solute Methanol =  */
/*   { */
/*     { "C"},   */
/*     { {-2.0,0.0, 0.0} */
/*     }, */
/*     {3.6705},        */
/*     {0.0800},        */
/*     {-0.040},         */
/*     1 */
/*   }; */






void RecomputeInitialSoluteData_Methanol(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Methanol.\n");
  RecomputeInitialSoluteData_II(BHD, &Methanol, damp, damp_LJ);
}



/*********************************/
/* Hexane */
/*********************************/
static Solute Hexane =
  {
    { "C","C","C","C","C","C","H","H","H","H","H","H","H","H","H","H","H","H","H","H"},
    { {1.709,  -2.812,   0.000},
      {1.684,  -1.278,   0.000},
      {0.245,  -0.753,   0.000},
      {0.241,   0.779,   0.000},
      {-1.198,   1.304,   0.000},
      {-1.206,   2.834,   0.000},
      {2.236,  -3.164,   0.887},
      {2.232,  -3.164,  -0.890},
      {0.691,  -3.204,   0.003},
      {2.202,  -0.914,  -0.888},
      {2.201,  -0.914,   0.890},
      {-0.273,  -1.115,   0.889},
      {-0.272,  -1.115,  -0.890},
      {0.757,   1.142,  -0.890},
      {0.757,   1.141,   0.890},
      {-1.716,   0.944,   0.890},
      {-1.716,   0.944,  -0.890},
      {-0.696,   3.204,  -0.890},
      {-0.696,   3.204,   0.890},
      {-2.236,   3.190,   0.000} },
    {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 2.5, 2.5, 2.5, 2.5,2.5, 2.5,2.5, 2.5,2.5, 2.5,2.5, 2.5,2.5, 2.5},
    {0.066, 0.066, 0.066, 0.066, 0.066, 0.066, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03},
    {-0.180, -0.120, -0.120, -0.120, -0.120, -0.180, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06 },
    20
  };

void RecomputeInitialSoluteData_Hexane(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Hexane.\n");
  RecomputeInitialSoluteData_II(BHD, &Hexane, damp, damp_LJ);
}

/* BUTANOIC ACID */
/* H1 sigma and epsilon adopted */
static Solute ButanoicAcid =
  {
    { "C1","O1","O2","C2","C3","C4", "OH", "H2", "H3", "H4", "H5", "H6", "H7", "H8"},
    {
      {1.422,  -0.017,   0.000},
      {1.422,   1.353,   0.000},
      {2.643,  -0.722,   0.000},
      {0.100,  -0.780,   0.000},
      {-1.060,   0.212,   0.000},
      {-2.381,  -0.551,   0.000},
      { 3.210,  -0.461,   0.882},
      { 0.043,  -1.407,   0.890},
      { 0.043,  -1.407,  -0.890},
      {-1.002,   0.838,  -0.890},
      {-1.002,   0.838,   0.890},
      {-2.439,  -1.178,   0.890},
      {-2.439,  -1.178,  -0.890},
      {-3.210,   0.157,   0.000} },
    {3.750, 2.960, 3.000, 3.500, 3.500, 3.500, 3.400, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500, 2.500},
    {0.105, 0.210, 0.170, 0.066, 0.066, 0.066, 0.046, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030},
    {0.520,-0.440,-0.530,-0.120,-0.120,-0.180, 0.450, 0.060, 0.060, 0.060, 0.060,  0.060, 0.060, 0.060},
    14
  };

void RecomputeInitialSoluteData_ButanoicAcid(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Butanoic Acid.\n");
  RecomputeInitialSoluteData_II(BHD, &ButanoicAcid, damp, damp_LJ);
}


/******************************/
/* Create initial solute data */
/******************************/
// XXX:  see (5.106) and (5.08) in the thesis
// XXX:  return BHD->gH_ini and BHD->gO_ini ( beta*(VM_LJ + VM_coulomb_short) )
// XXX:  and BHD->ucH, BHD->ucO ( beta * VM_coulomb_long ), but is beta missing here ??
static void RecomputeInitialSoluteData_II(BGY3dH2OData BHD, const Solute *S, real damp, real damp_LJ)
{
  DA da;
  PData PD;
  PetscScalar ***gHini_vec, ***gOini_vec;
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]);
  real r[3], r_s, h[3], interval[2], beta, L; //, fac;
  int x[3], n[3], i[3], N[3], k;
  real *p_O, *p_H;

  PD = BHD->PD;
  da = BHD->da;

  PetscPrintf(PETSC_COMM_WORLD,"Recomputing solute data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  p_O = (real*) malloc(3*sizeof(real));
  p_H = (real*) malloc(3*sizeof(real));

  VecSet(BHD->gH_ini, 0.0);
  VecSet(BHD->gO_ini, 0.0);
  VecSet(BHD->gHO_ini, 0.0);
  VecSet(BHD->ucH, 0.0);
  VecSet(BHD->ucO, 0.0);
  VecSet(BHD->ucHO, 0.0);
  FOR_DIM
    {
      VecSet(BHD->fH_l[dim],0.0);
      VecSet(BHD->fO_l[dim],0.0);
      VecSet(BHD->fHO_l[dim],0.0);
    }

  DAVecGetArray(da, BHD->gH_ini, &gHini_vec);
  DAVecGetArray(da, BHD->gO_ini, &gOini_vec);

  FOR_DIM
    {
      DAVecGetArray(da, BHD->fH_l[dim], &(fHl_vec[dim]));
      DAVecGetArray(da, BHD->fO_l[dim], &(fOl_vec[dim]));
    }

  /* loop over solute atoms */
  for(k=0; k<S->max_atoms; k++)
    {
      /* set parameters */
      p_O[0] = sqrt( eO * S->epsilon[k]);
      p_O[1] = 0.5*( sO + S->sigma[k]);
      p_O[2] =  qO * S->q[k];
      /* Empirical correction: add sigma 1.0 to small H */
      //if( !strcmp(S->names[k], "OH") )// && p_O[2]<0)
      //	p_O[1] += 1.0;

      p_H[0] = sqrt( eH * S->epsilon[k]);
      p_H[1] = 0.5*( sH + S->sigma[k]);
      p_H[2] =  qH * S->q[k];
      /* Empirical correction: add sigma 1.0 to small H */
      //if( !strcmp(S->names[k], "OH") )//&& p_H[2]<0)
      //	p_H[1] += 1.0;


      /* loop over local portion of grid */
      for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
	for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
	  for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	    {
	      /* set force vectors */
	      FOR_DIM
		r[dim] = i[dim]*h[dim]+interval[0] - S->x[k][dim];

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	      /* Lennard-Jones */
	      gHini_vec[i[2]][i[1]][i[0]] +=
		damp_LJ * beta* Lennard_Jones( r_s, p_H[0], p_H[1]);

	      gOini_vec[i[2]][i[1]][i[0]] +=
		damp_LJ * beta* Lennard_Jones( r_s, p_O[0], p_O[1]);

	      /* Coulomb short */
	      gHini_vec[i[2]][i[1]][i[0]] +=
		damp*beta* Coulomb_short( r_s, p_H[2]);

	      gOini_vec[i[2]][i[1]][i[0]] +=
		damp*beta* Coulomb_short( r_s, p_O[2]);
	    }

      ComputeSoluteDatafromCoulombII(BHD, BHD->v[0], S->x[k],  p_O[2], damp);
      VecAXPY(BHD->ucO, 1.0, BHD->v[0]);
      ComputeSoluteDatafromCoulombII(BHD, BHD->v[0], S->x[k],  p_H[2], damp);
      VecAXPY(BHD->ucH, 1.0, BHD->v[0]);

    }

  DAVecRestoreArray(da, BHD->gH_ini, &gHini_vec);
  DAVecRestoreArray(da, BHD->gO_ini, &gOini_vec);

  free(p_O);
  free(p_H);
}


#if 0
static void ComputeSoluteDatafromCoulomb(BGY3dH2OData BHD, Vec uc, const real x0[3], real q2, real damp)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], ic[3], N[3], index;
  real r[3], r_s, h[3], interval[2], k, fac, L, sign, fac2, L2;
  fftw_complex *fft_data;


  PD = BHD->PD;
  da = BHD->da;
  fft_data = BHD->g_fft;

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  L2=0.5*L;



  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));



  index=0;
   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0];


	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );


	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];

	    }
	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      fft_data[index].re = 0;
	      fft_data[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]))/SQR(L);
	      fac = EPSILON0INV/M_PI/k;
	      fac2= 2.*M_PI/L;
	      sign = cos(fac2*ic[0]*(x0[0]-L2))*cos(fac2*ic[1]*(x0[1]-L2))
		*cos(fac2*ic[2]*(x0[2]-L2));
	      //sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);
	      /* potential */
	      fft_data[index].re =   damp * q2 * sign *fac * exp(-k*SQR(M_PI)/SQR(G));
	      fft_data[index].im = 0;

	    }
	  index++;

	}

  /* FFT potential */
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, uc, fft_data,
			 BHD->fft_scratch, x, n, 0);
  VecScale(uc, 1./L/L/L);




/*  VecView(uc,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}
#endif


static void ComputeSoluteDatafromCoulombII(BGY3dH2OData BHD, Vec uc, const real x0[3], real q2, real damp)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], ic[3], N[3], index;
  real r[3], r_s, h[3], interval[2], k, fac, L, sign, h3;
  fftw_complex *fft_data, *dg_fft;
  PetscScalar ***v_vec;

  PD = BHD->PD;
  da = BHD->da;
  fft_data = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];
  h3 = h[0]*h[1]*h[2];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];




  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, uc, &v_vec);
  fac = pow(G/sqrt(M_PI),3.0);

   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0]-x0[dim];

          // XXX:  Gaussian distribution
	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  v_vec[i[2]][i[1]][i[0]] = fac*exp(-SQR(r_s*G));



	}
  DAVecRestoreArray(da, uc, &v_vec);
  // XXX:  convert to fft_complex
  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, uc, fft_data,
 			 BHD->fft_scratch, x, n, 0);


    index=0;
   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0];


	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );


	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];

	    }
	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft[index].re = 0;
	      dg_fft[index].im = 0;
	    }
	  else
          // XXX:  Fourier component of V_coulomb_long ??
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]))/SQR(L);
	      fac = h3*EPSILON0INV/M_PI/k;
	      sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);

	      dg_fft[index].re = damp * q2 *fac * fft_data[index].re;
	      dg_fft[index].im = damp * q2 *fac * fft_data[index].im;


	    }
	  index++;

	}

  /* FFT */
  // fft_complex to Vec
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, uc, dg_fft,
			 BHD->fft_scratch, x, n, 0);
  VecScale(uc, 1./L/L/L);




/*  VecView(uc,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

static void RecomputeInitialSoluteData_QM(BGY3dH2OData BHD, const Solute *S, real damp, real damp_LJ)
{
    PetscPrintf(PETSC_COMM_WORLD,"Recomputing solute(QM) data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);

    /*
     * Calculate LJ potential for all solvent sites.
     *
     * Beta is  the (inverse) temperature. For  historical reasons the
     * solute  field acting on  solvent sites  is defined  having this
     * factor.
     */
    real factor = damp_LJ * BHD->PD->beta;

    /*
     * Fill LJ-interaction of  H and O sites with  the solute into the
     * respective arrays.
     *
     * FIXME: LJ-parameters,  (eH, sH) and (eO,  sO), for H  and O are
     * #defined at some obscure place:
     *
     * We  supply ljc()  as a  callback function  that is  supposed to
     * compute the interaction  of a charged LJ solvent  site with the
     * solute.
     *
     * At  this place  the (short  range) Coulomb  interaction  of the
     * solvent  site  with  the  solute was  deliberately  omitted  by
     * specifying zero charge of the solvent site:
     */
    field (BHD->da, BHD->PD, S, eH, sH, 0.0, factor, ljc, BHD->gH_ini);
    field (BHD->da, BHD->PD, S, eO, sO, 0.0, factor, ljc, BHD->gO_ini);

    /*
     * Compute  the  charge  density  of  the  solute.   The  callback
     * function rho()  sums charge  distribution for each  solute site
     * and  does not use  (epsilon, sigma,  charge) parameters  of the
     * solvent site,  so that we  provide -1.0 for them.   The overall
     * factor is 1.0 (idependent of the solvent charge):
     */

    Vec rho_solute; /* Vector for solute charge density */

    DACreateGlobalVector (BHD->da, &rho_solute);

    field (BHD->da, BHD->PD, S, -1.0, -1.0, -1.0, 1.0, rho, rho_solute);

    /*
     * This solves  the Poisson equation and  puts resulting potential
     * into a pre-allocated (?) vector BHD->v[0].
     */
    ComputeSoluteDatafromCoulomb_QM (BHD, BHD->v[0], rho_solute, 1.0, damp);

    VecDestroy (rho_solute);

    VecSet (BHD->ucH, 0.0);
    VecAXPY (BHD->ucH, qH, BHD->v[0]);

    VecSet (BHD->ucO, 0.0);
    VecAXPY (BHD->ucO, qO, BHD->v[0]);
}

#if 0
/*
 * For debug purposes:
 */
static void dump (BGY3dH2OData BHD)
{
    PetscViewer viewer;

    /* print  to check */
    PetscPrintf(PETSC_COMM_WORLD, "Writing binarys \n");
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "LJH.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->gH_ini, viewer);
    PetscViewerDestroy(viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "LJO.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->gO_ini, viewer);
    PetscViewerDestroy(viewer);
    /* PetscViewerBinaryOpen(PETSC_COMM_WORLD, "gs.bin", FILE_MODE_WRITE, &viewer); */
    /* VecView(sumgs, viewer); */
    PetscViewerDestroy(viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "ucH.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->ucH, viewer);
    PetscViewerDestroy(viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "ucO.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->ucO, viewer);
    PetscViewerDestroy(viewer);
}
#endif

/*
 * Calculate a  real field "f"  for the solvent site  characterized by
 * (epsilon, sigma)  in the presence of  the solute S  with an overall
 * factor "fact" at every point (x, y, z) of the local grid.
 *
 * The function f(x, y, z, eps, sig,  chg, S) can be ljc() or rho() as
 * two examples.
 *
 * Vector "v" is the intent(out) argument.
 */
static void field (DA da, const ProblemData *PD,
                   const Solute *S,
                   real epsilon, real sigma, real charge, real fact,
                   real (*f)(real x, real y, real z,
                             real eps, real sig, real chg,
                             const Solute *S),
                   Vec v)
{
    PetscScalar ***vec;
    real h[3];
    int i0, j0, k0;
    int ni, nj, nk;

    /*
     * FIXME: do we really assume that intervals for x-, y- and z- are
     * the same? This  basically means the corner of  the unit cell is
     * at (offset, offset, offset):
     */
    real offset = PD->interval[0];

    FOR_DIM
        h[dim] = PD->h[dim];

    /* Get local portion of the grid */
    DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

    DAVecGetArray (da, v, &vec);

    /* loop over local portion of grid */
    for (int k = k0; k < k0 + nk; k++) {

        real z = k * h[2] + offset;

        for (int j = j0; j < j0 + nj; j++) {

            real y = j * h[1] + offset;

            for (int i = i0; i < i0 + ni; i++) {

                real x = i * h[0] + offset;

                /*
                 * Compute the field f at (x, y, z) <-> (i, j, k) e.g.
                 * by summing (LJ) contributions from all solute sites
                 * at that grid point:
                 */
                vec[k][j][i] = fact * f (x, y, z, epsilon, sigma, charge, S);
            }
        }
    }
    DAVecRestoreArray (da, v, &vec);
}

/*
 * Interaction of a charged LJ site (epsilon, sigma, charge) at (x, y,
 * z) with the solute S:
 */
static real ljc (real x, real y, real z,
                 real epsilon, real sigma, real charge,
                 const Solute *S)
{
    /* Sum force field contribution from all solute sites: */
    real field = 0.0;

    for (int site = 0; site < S->max_atoms; site++) {

        /* Interaction parameters for a pair of LJ sites: */
        real e2 = sqrt (epsilon * S->epsilon[site]);
        real s2 = 0.5 * (sigma + S->sigma[site]);

        /* Distance from a grid point to this site: */
        real r_s = sqrt (SQR(x - S->x[site][0]) +
                         SQR(y - S->x[site][1]) +
                         SQR(z - S->x[site][2]));

        /* 1. Lennard-Jones */
        field += Lennard_Jones (r_s, e2, s2);

        /* 2. Coulomb,  short range part.  For  historical reasons the
           overall scaling factor, the  product of solvent- and solute
           site charges, is handled by the function itself: */
        field += Coulomb_short (r_s, charge * S->q[site]);
    }

    return field;
}

/*
 * Charge  density  of  the solute  S  at  (x,  y, z).   Solvent  site
 * parameters (epsilon, sigma, charge)  are here to keep the interface
 * of rho() the same as that of ljc().
 *
 * Each gaussian is evaluated as:
 *
 *   rho(r) = q * [ G / sqrt(pi)]^3 * exp[-G^2 * (r - x0)^2]
 */
static real rho (real x, real y, real z,
                 real epsilon, real sigma, real charge, /* all three unused */
                 const Solute *S)
{
    /* G  is predefind  in bgy3d_SolventParameters.h  FIXME:  make the
       gaussian width a property of  the (solute) site in the same way
       as the charge of the site. */
    real prefac = pow(G / sqrt(M_PI), 3.0);

    /* Sum Gaussian contributions from all solute sites: */
    real field = 0.0;

    for (int site = 0; site < S->max_atoms; site++) {

        /* Square of the distance from a grid point to this site: */
        real r2 = (SQR(x - S->x[site][0]) +
                   SQR(y - S->x[site][1]) +
                   SQR(z - S->x[site][2]));

        /* Gaussian  distribution, note  that G  is not  a  width, but
           rather an inverse of it: */
        field += prefac * S->q[site] * exp(- G * G * r2);
    }

    return field;
}

// Solve Poisson Equation in Fourier space and get elestrostatic potential by inverse FFT
void ComputeSoluteDatafromCoulomb_QM(BGY3dH2OData BHD, Vec uc, Vec gs, real q, real damp)
{
    int x[3], n[3], i[3], ic[3], N[3], index;
    real h[3], interval[2], k, fac, L, h3; /* , sign; */
    fftw_complex *fft_gs, *fft_uc;

    interval[0] = BHD->PD->interval[0];
    interval[1] = BHD->PD->interval[1];
    L = interval[1] - interval[0];
    FOR_DIM
        h[dim] = BHD->PD->h[dim];
    FOR_DIM
        N[dim] = BHD->PD->N[dim];
    h3 = h[0] * h[1] * h[2];
    fft_gs = BHD->g_fft;
    fft_uc = BHD->gfg2_fft;
    VecSet(uc, 0.0);
    /* Get local portion of the grid */
    DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
    /* Get fft of gs
      gs(i, j, k) => fft_gs(ki, kj, kk)*/
    ComputeFFTfromVec_fftw(BHD->da, BHD->fft_plan_fw, gs, fft_gs, BHD->fft_scratch, x, n, 0);

    /* Solving Poisson Equation(SI) with FFT and IFFT:
      -1.0 * LAPLACIAN(V(x, y, z)) = 1/epsilon0 * rho(x, y, z)
    because: x = i * h, y = i * h, z = k * h, h = L / n
    =>-1.0 * n^2/L^2 * LAPLACIAN(uc(i, j, k)) = 1/epsilon0 * gs(i, j, k)
    FFT(see FFTW manual "What FFTW Really Computes"):
    => 4*pi^2 * n^2/L^2 * (kx^2 + ky^2 + kz^2)/n^2 * fft_uc(kx, ky, kz) = 1/epsilon0 * fft_gs(kx, ky, kz)
    => fft_uc(kx, ky, kz) = 1 / [(4 * pi * epsilon0) * pi * (kx^2 + ky^2 + kz^2)/L^2] * fft_gs(kx, ky, kz)
    IFFT(see FFTW manual "What FFTW Really Computes"):
    because: IFFT(fft_uc(kx, ky, kz)) = n^3 * uc(i, j, k)
    => uc(i, j, k) = h^3/L^3 * IFFT(fft_uc(kx, ky, kz)) */
    index = 0;
    for(i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    {
        for(i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
        {
            for(i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
            {

                FOR_DIM
                {
                    if( i[dim] <= N[dim] / 2)
                        ic[dim] = i[dim];
                    else
                        ic[dim] = i[dim] - N[dim];
                }

                if(ic[0] == 0 && ic[1] == 0 && ic[2] == 0)
                {
                    fft_uc[index].re = 0;
                    fft_uc[index].im = 0;
                }
                else
                {
                    k = (SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0])) / SQR(L);
                    // EPSILON0INV = 1 / 4 * pi * epsilon0
                    fac = h3 * EPSILON0INV / M_PI / k;
                    // sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);
                    fft_uc[index].re = damp * q * fac * fft_gs[index].re;
                    fft_uc[index].im = damp * q * fac * fft_gs[index].im;
                }
                index++;
            }
        }
    }
    /* IFFT(fft_uc(kx, ky, kz)) / L^3 */
    ComputeVecfromFFT_fftw(BHD->da, BHD->fft_plan_bw, uc, fft_uc, BHD->fft_scratch, x, n, 0);
    VecScale(uc, 1./L/L/L);
}
