/*==========================================================*/
/*  $Id: bgy3dtest.c,v 1.10 2007-03-29 07:44:34 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fft.h"
#include "bgy3ddiv.h"
#include "bgy3dtest.h"
#include "bgy3dfourier.h"

static void ComputeError(Vec gmax, BGY3dFourierData BDDmax, int Nmax, Vec g, BGY3dFourierData BDD, int N);

#define NMAX 128

Vec BGY3dDiv_solve_FourierTest(ProblemData *PD, Vec g_ini, int vdim)
{
  Vec g[5];
  BGY3dFourierData BDD[5];
  int i, index;


  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Discretization test for BGY3dDiv equation with Fourier ansatz...\n");


  index=0;
  for(i=NMAX; i>=32 ; i/=2)
    {
      FOR_DIM
	PD->N[dim]=i;
      FOR_DIM
	PD->h[dim] = (PD->interval[1]-PD->interval[0])/i;

      PetscPrintf(PETSC_COMM_WORLD, "===============================\n");
      PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
      PetscPrintf(PETSC_COMM_WORLD, "N = %d\n", PD->N[0]);
      g[index]=BGY3dDiv_solve_Fourier(PD, PETSC_NULL, 0);
      BDD[index++] = BGY3dFourierData_kirk_malloc(PD);
    }
  FOR_DIM
    BDD[0]->PD->N[dim]=NMAX;


  /***************************/
  /* Factor e^-v at highest level */
  /***************************/
  //VecCopy(g[0], BDD[0]->v[0]);
  //Compute_g(BDD[0], g[0], BDD[0]->g_ini, BDD[0]->v[0]);
  /***************************/
  ExtractAxis(BDD[0], g[0], 0);
  index=1;
  for(i=NMAX/2; i>=32 ; i/=2)
    {
      ComputeError(g[0], BDD[0], NMAX, g[index], BDD[index], i);
      index++;
    }
  index=0;
  for(i=NMAX; i>=32 ; i/=2)
    {
      BGY3dFourierData_free(BDD[index]);
      VecDestroy(g[index]);
      index++;
    }



  return PETSC_NULL;
}


static void ComputeError(Vec gmax, BGY3dFourierData BDDmax, int Nmax, Vec g, BGY3dFourierData BDD, int N)
{
  DA damax, da;
  Vec dgint, gint;
  PetscScalar ***gint_vec, ***g_vec, norm2, norminf;
  int i[3], n[3], x[3], il[3], ili[3], gi[3];
  real h, v[3];


  PetscPrintf(PETSC_COMM_WORLD,"Interpolating from %d to %d...\n", N, Nmax);

  FOR_DIM
    BDD->PD->N[dim]=Nmax;
  damax=BDDmax->da;
  da = BDD->da;
  h=Nmax/N;

  VecDuplicate(gmax, &gint);
  VecDuplicate(gmax, &dgint);


  DAGetCorners(damax, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  VecSet(gint,0.0);
  DAVecGetArray(damax, gint, &gint_vec);
  DAVecGetArray(da, g, &g_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  FOR_DIM
	    il[dim] = i[dim]/h;

	  for(ili[2]=0; ili[2]<2; ili[2]++)
	    {
	      v[2] = 1.-abs(i[2]-(il[2]+ili[2])*h)/h;

	      for(ili[1]=0; ili[1]<2; ili[1]++)
		{
		  v[1]= 1.-abs(i[1]-(il[1]+ili[1])*h)/h;

		  for(ili[0]=0; ili[0]<2; ili[0]++)
		    {
		      v[0]= 1.-abs(i[0]-(il[0]+ili[0])*h)/h;


		      FOR_DIM
			{
			  gi[dim] = il[dim]+ili[dim];
			  if(gi[dim]>= N )
			    gi[dim]=0;
			}
		      gint_vec[i[2]][i[1]][i[0]] +=
			v[2]*v[1]*v[0]*g_vec[gi[2]][gi[1]][gi[0]];
		    }
		}
	    }
	}
  DAVecRestoreArray(damax, gint, &gint_vec);
  DAVecRestoreArray(da, g, &g_vec);

  /***************************/
  /* Factor e^-v at highest level */
  /***************************/
  //VecCopy(gint, dgint);
  //Compute_g(BDDmax, gint, BDDmax->g_ini, dgint);
  /***************************/

/*   VecView(gint,PETSC_VIEWER_STDERR_WORLD);   */
  ExtractAxis(BDDmax, gint, 0);

  VecAXPY(gint, -1.0, gmax);
  VecNorm(gint, NORM_2, &norm2);
  VecNorm(gint, NORM_INFINITY, &norminf);

  PetscPrintf(PETSC_COMM_WORLD,"L2= %e \t Linfty= %e \n", norm2/pow(Nmax,3.0), norminf);

  VecDestroy(dgint);
  VecDestroy(gint);
}


Vec BGY3dDiv_test(ProblemData *PD, Vec g_ini, int vdim)
{
  BGY3dDivData BDD;
  Vec g,  f, rhs;
  Mat SM;
  PetscScalar sigma_g=1.0, sigma_K=1.0, f_norm_l2, f_norm_max;
  int N3;
  PetscTruth flg;




  PetscPrintf(PETSC_COMM_WORLD, "Testing BGY3dDiv model...\n");

  flg = bgy3d_getopt_test ("-seq");
  BDD = BGY3dDivData_malloc(PD,flg);
  if(flg)
    DAGetMatrix( BDD->da, MATSEQAIJ, &SM);
  else
    DAGetMatrix( BDD->da, MATMPIAIJ, &SM);

  DACreateGlobalVector(BDD->da, &g);
  DACreateGlobalVector(BDD->da, &rhs);
  DACreateGlobalVector(BDD->da, &f);


  /* read BGY3dDiv specific things from command line */
  /* Mixing parameter */
  bgy3d_getopt_real ("-sigma_g", &sigma_g);
  bgy3d_getopt_real ("-sigma_K", &sigma_K);


  InitializeTestData(BDD, g, sigma_g, sigma_K);
  ComputeIntegralPart(BDD, g, f);

  /* Assemble matrix to solve */
  AssembleSystemMatrix(BDD, SM, f);
  AssembleSystemMatrix_part2b(BDD, SM);

  /* set right hand site */
  ComputeRHStest(BDD, g, rhs, sigma_g, sigma_K);


  MatMult(SM, g, f);
  VecAXPY(f,-1.0, rhs);
  //VecView(rhs,PETSC_VIEWER_STDERR_WORLD);
  //VecWAXPY(f, -1.0, BDD->i[0], rhs);

  N3 = PD->N[0]*PD->N[1]*PD->N[2];
  VecNorm(f, NORM_2, &f_norm_l2);
  VecNorm(f, NORM_INFINITY, &f_norm_max);
  PetscPrintf(PETSC_COMM_WORLD,"L2-norm of error: %e\n",f_norm_l2/N3);
  PetscPrintf(PETSC_COMM_WORLD,"Max-norm of error: %e\n",f_norm_max);

  VecDestroy(rhs);
  VecDestroy(g);
  MatDestroy(SM);

  BGY3dDivData_free(BDD);

  return f;

}

Vec BGY3dDivFourier_test(ProblemData *PD, Vec g_ini, int vdim)
{
  BGY3dDivData BDD;
  Vec g,  f, rhs;
  PetscScalar sigma_g=1.0, sigma_K=1.0, f_norm_l2, f_norm_max;
  int N3;
  PetscTruth flg;




  PetscPrintf(PETSC_COMM_WORLD, "Testing BGY3dDivFourier model...\n");

  BDD = BGY3dDivData_malloc(PD,flg);

  DACreateGlobalVector(BDD->da, &g);
  DACreateGlobalVector(BDD->da, &rhs);
  DACreateGlobalVector(BDD->da, &f);


  /* read BGY3dDiv specific things from command line */
  /* Mixing parameter */
  bgy3d_getopt_real ("-sigma_g", &sigma_g);
  bgy3d_getopt_real ("-sigma_K", &sigma_K);


  InitializeTestData(BDD, g, sigma_g, sigma_K);
  VecShift(g, 1.0); /* g=g-1 in Compute_dg ! */
  //Compute_dg(BDD,  g,  f);

  /* set right hand site */
  //ComputeRHStestFourier(BDD, g, rhs, sigma_g, sigma_K);


  VecAXPY(f,-1.0, rhs);
  //VecView(rhs,PETSC_VIEWER_STDERR_WORLD);
  //VecWAXPY(f, -1.0, BDD->i[0], rhs);

  N3 = PD->N[0]*PD->N[1]*PD->N[2];
  VecNorm(f, NORM_2, &f_norm_l2);
  VecNorm(f, NORM_INFINITY, &f_norm_max);
  PetscPrintf(PETSC_COMM_WORLD,"L2-norm of error: %e\n",f_norm_l2/N3);
  PetscPrintf(PETSC_COMM_WORLD,"Max-norm of error: %e\n",f_norm_max);

  VecDestroy(rhs);
  VecDestroy(g);

  BGY3dDivData_free(BDD);

  return f;

}




void InitializeTestData(BGY3dDivData BDD, Vec g, real sigma_g, real sigma_K)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], N[3];
  PetscScalar ***g_vec, ***(k_vec[3]);
  real h[3], r[3], r_s, facg, facK, L;

  da = BDD->da;
  PD = BDD->PD;

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];
  L = PD->interval[1]-PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, g, &(g_vec));
  FOR_DIM
    DAVecGetArray(da, BDD->v[dim], &(k_vec[dim]));

  facg = 1./pow(2.*M_PI*SQR(sigma_g),1.5);
  facK = 1./(sigma_K*sqrt(2.*M_PI));
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+PD->interval[0];
	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );


	  g_vec[i[2]][i[1]][i[0]] = facg * exp(-0.5/SQR(sigma_g)*SQR(r_s));


	  FOR_DIM
	    {
	      r[dim] = i[dim]*h[dim];
	      if( i[dim]>=N[dim]/2)
		r[dim] -= L;
	    }

	  FOR_DIM
	    k_vec[dim][i[2]][i[1]][i[0]] =
	    facK * exp(-0.5/SQR(sigma_K)*SQR(r[dim]));

	}
  DAVecRestoreArray(da, g, &(g_vec));
  FOR_DIM
    DAVecRestoreArray(da, BDD->v[dim], &(k_vec[dim]));

  FOR_DIM
    BDD->fg2_fft[dim] = ComputeFFTfromVec(da, BDD->fft_plan,
					  BDD->v[dim],
					  BDD->fg2_fft[dim]);


}


void ComputeRHStest(BGY3dDivData BDD, Vec g, Vec rhs, real sigma_g,
		    real sigma_K)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], N[3];
  PetscScalar ***g_vec, ***rhs_vec;
  real h[3], r[3], r_s, facg, facconv, L, s2, ss2;

  da = BDD->da;
  PD = BDD->PD;

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];
  L = PD->interval[1]-PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, g, &(g_vec));
  DAVecGetArray(da, rhs, &(rhs_vec));


  facg = 1./(sigma_g*sqrt(2.*M_PI));
  s2=1./SQR(sigma_g);
  ss2 = 1./(SQR(sigma_g)+SQR(sigma_K));
  facconv = PD->beta*PD->rho/(sqrt(2.*M_PI*(SQR(sigma_g)+SQR(sigma_K)))) *
    (s2+ss2);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+PD->interval[0];
	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );



	  rhs_vec[i[2]][i[1]][i[0]] = g_vec[i[2]][i[1]][i[0]] *
	    ( SQR(s2)*SQR(r_s) - s2*3.0 -  facconv *
	      ( exp(-SQR(r[0])/(2.*(SQR(sigma_g)+SQR(sigma_K)))) * r[0] +
		exp(-SQR(r[1])/(2.*(SQR(sigma_g)+SQR(sigma_K)))) * r[1] +
		exp(-SQR(r[2])/(2.*(SQR(sigma_g)+SQR(sigma_K)))) * r[2] ));

/* 	   rhs_vec[i[2]][i[1]][i[0]] =  */
/* 	     PD->beta*PD->rho/(sqrt(2.*M_PI*(SQR(sigma_g)+SQR(sigma_K)))) * */
/* 	     exp(-SQR(r[0])/(2.*(SQR(sigma_g)+SQR(sigma_K)))); */

	}
  DAVecRestoreArray(da, g, &(g_vec));
  DAVecRestoreArray(da, rhs, &(rhs_vec));

}

void ComputeRHStestFourier(BGY3dDivData BDD, Vec g, Vec rhs, real sigma_g,
			   real sigma_K)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], N[3];
  PetscScalar ***g_vec, ***rhs_vec;
  real h[3], r[3], r_s, facg, facconv, L, s2, ss2;

  da = BDD->da;
  PD = BDD->PD;

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];
  L = PD->interval[1]-PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, g, &(g_vec));
  DAVecGetArray(da, rhs, &(rhs_vec));


  facg = 1./(sigma_g*sqrt(2.*M_PI));
  s2=1./SQR(sigma_g);
  ss2 = 1./(SQR(sigma_g)+SQR(sigma_K));
  facconv = PD->beta*PD->rho/(sqrt(2.*M_PI*(SQR(sigma_g)+SQR(sigma_K)))) *
    (s2+ss2);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+PD->interval[0];
	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );



	  rhs_vec[i[2]][i[1]][i[0]] = g_vec[i[2]][i[1]][i[0]] *
	    ( SQR(s2)*SQR(r_s) - s2*3.0 -  facconv *
	      ( exp(-SQR(r[0])/(2.*(SQR(sigma_g)+SQR(sigma_K)))) * r[0] +
		exp(-SQR(r[1])/(2.*(SQR(sigma_g)+SQR(sigma_K)))) * r[1] +
		exp(-SQR(r[2])/(2.*(SQR(sigma_g)+SQR(sigma_K)))) * r[2] ));



	}
  DAVecRestoreArray(da, g, &(g_vec));
  DAVecRestoreArray(da, rhs, &(rhs_vec));

}


void InitializeConvolutionData(BGY3dFourierData BDD, real sigma_g1, real sigma_g2,
			       Vec gg, Vec sol)
{
  ProblemData *PD;
  Vec g1, g2;
  DA da;
  int x[3], n[3], i[3], N[3], k;
  PetscScalar ***g1_vec, ***g2_vec, ***s_vec;
  real h[3], r[3], r_s, facg1, facg2, facs, L, max_k;
  FFT_DATA *(fg2_fft[3]), *g_fft;

  da = BDD->da;
  PD = BDD->PD;

  g1 = BDD->f[0];
  g2 = BDD->f[1];

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];
  L = PD->interval[1]-PD->interval[0];

  g_fft = BDD->g_fft;
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, g1, &(g1_vec));
  DAVecGetArray(da, g2, &(g2_vec));
  DAVecGetArray(da, sol, &(s_vec));

  facg1 = 1./pow(2.*M_PI*SQR(sigma_g1),1.5);
  facg2 = 1./pow(2.*M_PI*SQR(sigma_g2),1.5);
  facs  = 1./pow(2.*M_PI*(SQR(sigma_g1)+SQR(sigma_g2)),1.5);


  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+PD->interval[0];
	  r_s = ( SQR(r[0])+SQR(r[1])+SQR(r[2]) );


	  g1_vec[i[2]][i[1]][i[0]] = facg1 * exp(-0.5/SQR(sigma_g1)*r_s);
	  s_vec[i[2]][i[1]][i[0]]  = facs * exp(-0.5/(SQR(sigma_g1)+SQR(sigma_g2))*r_s);

	  FOR_DIM
	    {
	      r[dim] = i[dim]*h[dim];
	      if( i[dim]>=N[dim]/2)
		r[dim] -= L;
	    }
	  r_s =  SQR(r[0])+SQR(r[1])+SQR(r[2]) ;

	  g2_vec[i[2]][i[1]][i[0]] = facg2 * exp(-0.5/SQR(sigma_g2)*r_s);

	}
  DAVecRestoreArray(da, g1, &(g1_vec));
  DAVecRestoreArray(da, g2, &(g2_vec));
  DAVecRestoreArray(da, sol, &(s_vec));


  /* Compute convolution */
  ComputeFFTfromVec(da, BDD->fft_plan, g1, fg2_fft[0]);
  ComputeFFTfromVec(da, BDD->fft_plan, g2, fg2_fft[1]);

  max_k=n[0]*n[1]*n[2];
  for(k=0; k<max_k; k++)
    {
      g_fft[k].re = fg2_fft[0][k].re*fg2_fft[1][k].re
	       - fg2_fft[0][k].im*fg2_fft[1][k].im;

      g_fft[k].im = fg2_fft[0][k].re*fg2_fft[1][k].im
	       + fg2_fft[0][k].im*fg2_fft[1][k].re;
    }

  ComputeVecfromFFT(da, BDD->fft_plan, gg, g_fft);

  VecScale(gg, PD->h[0]*PD->h[1]*PD->h[2]/PD->N[0]/PD->N[1]/PD->N[2]);


}


Vec BGY3d_Convolution_Test(ProblemData *PD, Vec g_ini, int vdim)
{
  BGY3dFourierData BDD;
  Vec sol,  f, gg;
  PetscScalar sigma_g1=1.0, sigma_g2=1.0, f_norm_l2, f_norm_max;
  int N3;

  PetscPrintf(PETSC_COMM_WORLD, "Testing Convolution ...\n");

  BDD = BGY3dFourierData_kirk_malloc(PD);

  gg = BDD->f[2];
  sol = BDD->g_ini;
  f = BDD->v[2];



  /* read BGY3dDiv specific things from command line */
  /* Mixing parameter */
  bgy3d_getopt_real ("-sigma_g1", &sigma_g1);
  bgy3d_getopt_real ("-sigma_g2", &sigma_g2);
  PetscPrintf(PETSC_COMM_WORLD,"sigma_g1= %f\n", sigma_g1);
  PetscPrintf(PETSC_COMM_WORLD,"sigma_g2= %f\n", sigma_g2);

  InitializeConvolutionData(BDD, sigma_g1, sigma_g2, gg, sol);






  VecWAXPY(f,-1.0, sol, gg);

  /* VecView(gg,PETSC_VIEWER_STDERR_WORLD);   */

  N3 = PD->N[0]*PD->N[1]*PD->N[2];
  VecNorm(f, NORM_2, &f_norm_l2);
  VecNorm(f, NORM_INFINITY, &f_norm_max);
  PetscPrintf(PETSC_COMM_WORLD,"L2-norm of error: %e\n",f_norm_l2/N3);
  PetscPrintf(PETSC_COMM_WORLD,"Max-norm of error: %e\n",f_norm_max);



  BGY3dFourierData_free(BDD);

  return PETSC_NULL;

}
