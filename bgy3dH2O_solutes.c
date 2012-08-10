/*==========================================================*/
/*  $Id: bgy3dH2O_solutes.c,v 1.3 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"

#include "bgy3d_SolventParameters.h"



void ComputeSoluteDatafromCoulomb(BGY3dH2OData BHD, Vec uc, real x0[3], real q2,
				  real damp);
void ComputeSoluteDatafromCoulombII(BGY3dH2OData BHD, Vec uc, real x0[3], real q2,
				  real damp);

#define MAXATOM 20

typedef struct SoluteStruct
{
  char names[MAXATOM][5];       /* atom types */
  real x[MAXATOM][3];           /* the coordinates  */
  real sigma[MAXATOM];          /* sigma for LJ */
  real epsilon[MAXATOM];        /* epsilon for LJ */
  real q[MAXATOM];              /* charges */
  int max_atoms;
}Solute;

void ComputeSoluteDatafromLJ_QM(BGY3dH2OData BHD, Solute *S, real damp_LJ);
void ComputeSoluteDatafromCoulomb_QM(BGY3dH2OData BHD, Vec uc, Vec gs, real q, real damp);
void CreateGaussian(BGY3dH2OData BHD, Vec gs, real q, real widht, real *x0);
// void ConvertQMSolute(Solute *S, QM_Solute *QMS);
void RecomputeInitialSoluteData_QM(BGY3dH2OData BHD, Solute *S, real damp, real damp_LJ, real zpad);

void RecomputeInitialSoluteData_II(BGY3dH2OData BHD, Solute *S, real damp, real damp_LJ, real zpad);

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
  RecomputeInitialSoluteData_II(BHD, &Water, damp, damp_LJ, zpad);
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
  RecomputeInitialSoluteData_II(BHD, &CarbonDisulfide, damp, damp_LJ, zpad);
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
  RecomputeInitialSoluteData_II(BHD, &HydrogenChloride, damp, damp_LJ, zpad);
#endif
#ifdef QM
  RecomputeInitialSoluteData_QM(BHD, &HydrogenChloride, damp, damp_LJ, zpad);
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
  RecomputeInitialSoluteData_II(BHD, &Methanol, damp, damp_LJ, zpad);
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
  RecomputeInitialSoluteData_II(BHD, &Hexane, damp, damp_LJ, zpad);
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
  RecomputeInitialSoluteData_II(BHD, &ButanoicAcid, damp, damp_LJ, zpad);
}


/******************************/
/* Create initial solute data */
/******************************/
// XXX:  see (5.106) and (5.08) in the thesis
// XXX:  return BHD->gH_ini and BHD->gO_ini ( beta*(VM_LJ + VM_coulomb_short) )
// XXX:  and BHD->ucH, BHD->ucO ( beta * VM_coulomb_long ), but is beta missing here ??
void RecomputeInitialSoluteData_II(BGY3dH2OData BHD, Solute *S, real damp, real damp_LJ, real zpad)
{
  DA da;
  PData PD;
  PetscScalar ***gHini_vec, ***gOini_vec;
  // PetscScalar ***ucH_vec, ***ucO_vec;
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]);
  real r[3], r_s, h[3], interval[2], beta, L; //, fac;
  int x[3], n[3], i[3], dim, N[3], k;
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

/*   DAVecGetArray(da, BHD->ucH, &ucH_vec); */
/*   DAVecGetArray(da, BHD->ucO, &ucO_vec); */
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
		// damp_LJ * beta* Lennard_Jones( r_s, (void*)p_H);
		damp_LJ * beta* Lennard_Jones( r_s, p_H[0], p_H[1]);
	      gOini_vec[i[2]][i[1]][i[0]] +=
		// damp_LJ * beta* Lennard_Jones( r_s, (void*)p_O);
		damp_LJ * beta* Lennard_Jones( r_s, p_O[0], p_O[1]);

	      /* Coulomb short */
	      gHini_vec[i[2]][i[1]][i[0]] +=
		// damp*beta* Coulomb_short( r_s, (void*)p_H);
                // replace void pointer as real number
		damp*beta* Coulomb_short( r_s, p_H[2]);
	      gOini_vec[i[2]][i[1]][i[0]] +=
		// damp*beta* Coulomb_short( r_s, (void*)p_O);
                // replace void pointer as real number
		damp*beta* Coulomb_short( r_s, p_O[2]);


	      /* Coulomb long */
	      /* 	  gHini_vec[i[2]][i[1]][i[0]] +=  */
	      /* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsH); */
	      /* 	  gOini_vec[i[2]][i[1]][i[0]] +=  */
	      /* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */
	      /* 	  gHOini_vec[i[2]][i[1]][i[0]] +=  */
	      /* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */

	      /* Coulomb */
	      /* 	  gHini_vec[i[2]][i[1]][i[0]] +=  */
	      /* 	    damp * beta* Coulomb( r_s, BHD->LJ_paramsH); */
	      /* 	  gOini_vec[i[2]][i[1]][i[0]] +=  */
	      /* 	    damp * beta* Coulomb( r_s, BHD->LJ_paramsO); */
	      /* 	  gHOini_vec[i[2]][i[1]][i[0]] +=  */
	      /* 	    damp * beta* Coulomb( r_s, BHD->LJ_paramsHO); */

/* 	      ucH_vec[i[2]][i[1]][i[0]] +=  */
/* 		damp * beta* Coulomb_long( r_s, (void*)p_H); */
/* 	      ucO_vec[i[2]][i[1]][i[0]] +=  */
/* 		damp * beta* Coulomb_long( r_s, (void*)p_O); */



	    }

      ComputeSoluteDatafromCoulombII(BHD, BHD->v[0], S->x[k],  p_O[2], damp);
      VecAXPY(BHD->ucO, 1.0, BHD->v[0]);
      ComputeSoluteDatafromCoulombII(BHD, BHD->v[0], S->x[k],  p_H[2], damp);
      VecAXPY(BHD->ucH, 1.0, BHD->v[0]);

    }

  DAVecRestoreArray(da, BHD->gH_ini, &gHini_vec);
  DAVecRestoreArray(da, BHD->gO_ini, &gOini_vec);
/*   DAVecRestoreArray(da, BHD->ucH, &ucH_vec); */
/*   DAVecRestoreArray(da, BHD->ucO, &ucO_vec); */


/*   ComputeFFTSoluteII(BHD, BHD->ucH , BHD->ucHO, BHD->LJ_paramsHO, damp, zpad); */
/*   VecScale(BHD->ucH, beta); */
/*   VecAXPY( BHD->gH_ini, beta, BHD->ucHO); */

/*   ComputeFFTSoluteII(BHD, BHD->ucO , BHD->ucHO, BHD->LJ_paramsO,  damp, zpad); */
/*   VecAXPY( BHD->gO_ini, beta, BHD->ucHO); */
/*   VecScale(BHD->ucO, beta); */

  /* Shift uc's */
/*   VecSum(BHD->ucH, &fac); */
/*   VecShift(BHD->ucH, -fac/N[0]/N[1]/N[2]); */
/*   VecSum(BHD->ucO, &fac); */
/*   VecShift(BHD->ucO, -fac/N[0]/N[1]/N[2]); */

/*   VecView(BHD->ucH,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

  free(p_O);
  free(p_H);

}


void ComputeSoluteDatafromCoulomb(BGY3dH2OData BHD, Vec uc, real x0[3], real q2,
				  real damp)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], ic[3], N[3], dim, index;
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


void ComputeSoluteDatafromCoulombII(BGY3dH2OData BHD, Vec uc, real x0[3], real q2,
				    real damp)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], ic[3], N[3], dim, index;
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

void RecomputeInitialSoluteData_QM(BGY3dH2OData BHD, Solute *S, real damp, real damp_LJ, real zpad)
{
    DA da;
    Vec gs, sumgs; /* Vector for gaussian */
    int k;
    PetscViewer viewer;

    da = BHD->da;
    DACreateGlobalVector(da, &gs);
    DACreateGlobalVector(da, &sumgs);
    VecSet(gs, 0.0);
    VecSet(sumgs, 0.0);
    VecSet(BHD->ucH, 0.0);
    VecSet(BHD->ucO, 0.0);

    PetscPrintf(PETSC_COMM_WORLD,"Recomputing solute(QM) data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);
    // LJ potential energy
    ComputeSoluteDatafromLJ_QM(BHD, S, damp_LJ);

    // Create charge distribution for each atom center
    // then sum them up
    for(k=0; k<S->max_atoms; k++)
    {
        // G is predefind in bgy3d_SolventParameters.h
        CreateGaussian(BHD, gs, S->q[k], G, S->x[k]);
        VecAXPY(sumgs, 1.0, gs);
    }
    /*    ComputeSoluteDatafromCoulomb_QM(BHD, BHD->v[0], gs, qH, damp);
        VecAXPY(BHD->ucH, 1.0, BHD->v[0]);
        ComputeSoluteDatafromCoulomb_QM(BHD, BHD->v[0], gs, qO, damp);
        VecAXPY(BHD->ucO, 1.0, BHD->v[0]);
     } */
    ComputeSoluteDatafromCoulomb_QM(BHD, BHD->v[0], sumgs, qH, damp);
    VecAXPY(BHD->ucH, 1.0, BHD->v[0]);
    ComputeSoluteDatafromCoulomb_QM(BHD, BHD->v[0], sumgs, qO, damp);
    VecAXPY(BHD->ucO, 1.0, BHD->v[0]);
    /* print  to check */
    PetscPrintf(PETSC_COMM_WORLD, "Writing binarys \n");
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "LJH.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->gH_ini, viewer);
    PetscViewerDestroy(viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "LJO.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->gO_ini, viewer);
    PetscViewerDestroy(viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "gs.bin", FILE_MODE_WRITE, &viewer);
    VecView(sumgs, viewer);
    PetscViewerDestroy(viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "ucH.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->ucH, viewer);
    PetscViewerDestroy(viewer);
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "ucO.bin", FILE_MODE_WRITE, &viewer);
    VecView(BHD->ucO, viewer);
    PetscViewerDestroy(viewer);

    // exit(1);
    VecDestroy(gs);
    VecDestroy(sumgs);

}

// Calculate LJ potential
void ComputeSoluteDatafromLJ_QM(BGY3dH2OData BHD, Solute *S, real damp_LJ)
{
    PetscScalar ***gHini_vec, ***gOini_vec;
    real ffpara_H[2], ffpara_O[2];
    real r[3], r_s, interval[2], beta, h[3];
    int x[3], n[3], i[3], dim, site;

    interval[0] = BHD->PD->interval[0];
    FOR_DIM
        h[dim] = BHD->PD->h[dim];
    beta = BHD->PD->beta;
    /* Get local portion of the grid */
    DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

    VecSet(BHD->gH_ini, 0.0);
    VecSet(BHD->gO_ini, 0.0);
    DAVecGetArray(BHD->da, BHD->gH_ini, &gHini_vec);
    DAVecGetArray(BHD->da, BHD->gO_ini, &gOini_vec);

    /* loop over solute atoms, but only LJ interactions */
    for(site = 0; site < S->max_atoms; site++)
    {
        ffpara_O[0] = sqrt( eO * S->epsilon[site]);
        ffpara_O[1] = 0.5*( sO + S->sigma[site]);

        ffpara_H[0] = sqrt( eH * S->epsilon[site]);
        ffpara_H[1] = 0.5*( sH + S->sigma[site]);

        /* loop over local portion of grid */
        for(i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
        {
            for(i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
            {
                for(i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
                {
                    FOR_DIM
                        r[dim] = i[dim] * h[dim] + interval[0] - S->x[site][dim];

                    r_s = sqrt( SQR(r[0]) + SQR(r[1]) + SQR(r[2]));

                    /* Lennard-Jones */
                    gHini_vec[i[2]][i[1]][i[0]] +=
                         damp_LJ * beta * Lennard_Jones(r_s, ffpara_H[0], ffpara_H[1]);
                    gOini_vec[i[2]][i[1]][i[0]] +=
                         damp_LJ * beta * Lennard_Jones(r_s, ffpara_O[0], ffpara_O[1]);
                }
            }
        }
    }

    DAVecRestoreArray(BHD->da, BHD->gH_ini, &gHini_vec);
    DAVecRestoreArray(BHD->da, BHD->gO_ini, &gOini_vec);

}

// Create gaussian on 3d cartesian grid
// rho(r) = q * [ width / sqrt(pi)]^3 * exp[-width^2 * (r - x0)^2]
void CreateGaussian(BGY3dH2OData BHD, Vec gs, real q, real width, real *x0)
{
    PetscScalar ***gs_vec;
    real r[3], r_s, interval[2], h[3], prefac;
    int x[3], n[3], i[3], dim;

    interval[0] = BHD->PD->interval[0];
    FOR_DIM
        h[dim] = BHD->PD->h[dim];

    DAVecGetArray(BHD->da, gs, &gs_vec);

    /* Get local portion of the grid */
    DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

    prefac = pow(width / sqrt(M_PI), 3.0) * q;
    /* loop over local portion of grid */
    for(i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    {
        for(i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
        {
            for(i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
            {
                FOR_DIM
                    r[dim] = i[dim] * h[dim] + interval[0] - x0[dim];

                r_s =  SQR(r[0]) + SQR(r[1]) + SQR(r[2]);
                gs_vec[i[2]][i[1]][i[0]] = prefac * exp(-1.0 * width * width * r_s);
            }
        }
    }

    DAVecRestoreArray(BHD->da, gs, &gs_vec);
}

// Solve Poisson Equation in Fourier space and get elestrostatic potential by inverse FFT
void ComputeSoluteDatafromCoulomb_QM(BGY3dH2OData BHD, Vec uc, Vec gs, real q, real damp)
{
    int x[3], n[3], i[3], ic[3], N[3], dim, index;
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
