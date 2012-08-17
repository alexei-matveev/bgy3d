/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.15 2007-04-23 16:55:13 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"
#include "bgy3d-getopt.h"


#define rab 1.0



BGY3dDiatomicABData BGY3dDiatomicABData_Pair_malloc(ProblemData *PD)
{
  BGY3dDiatomicABData BDD;
  DA da;
  real interval[2], h[3], N[3], L, r[3], r_s, beta;
  int i[3], x[3], n[3];
  PetscScalar ***(fa_vec[3]),***(fb_vec[3]),***(fab_vec[3]);
  PetscScalar ***gaini_vec, ***gbini_vec, ***gabini_vec;
  int np;
  int local_nx, local_x_start, local_ny, local_y_start, total_local_size;
  PetscInt lx[1], ly[1], *lz;
  real epsilona, epsilonb, epsilonab;
  real sigmaa, sigmab, sigmaab;


  real eps[]={ 1.0, 1.0}, sig[]={ 0.5, 1.0};
  //real eps[]={ 0.1512, 0.046}, sig[]={ 3.1506, 0.4};


  BDD = (BGY3dDiatomicABData) malloc(sizeof(*BDD));


  // BDD->LJ_paramsa = (void* ) malloc(sizeof(real)*2);
  // ((real*)(BDD->LJ_paramsa))[0] = eps[0];   /* espilon */
  // ((real*)(BDD->LJ_paramsa))[1] = sig[0];   /* sigma   */
  BDD->LJ_paramsa[0] = eps[0];
  BDD->LJ_paramsa[1] = sig[0];
  epsilona = BDD->LJ_paramsa[0];
  sigmaa = BDD->LJ_paramsa[1];

  // BDD->LJ_paramsb = (void* ) malloc(sizeof(real)*2);
  // ((real*)(BDD->LJ_paramsb))[0] = eps[1];   /* espilon */
  // ((real*)(BDD->LJ_paramsb))[1] = sig[1];   /* sigma   */
  BDD->LJ_paramsb[0] = eps[1];
  BDD->LJ_paramsb[1] = sig[1];
  epsilonb = BDD->LJ_paramsb[0];
  sigmab = BDD->LJ_paramsb[1];

  // BDD->LJ_paramsab = (void* ) malloc(sizeof(real)*2);
  // ((real*)(BDD->LJ_paramsab))[0] = sqrt(eps[0]*eps[1]);   /* espilon */
  // ((real*)(BDD->LJ_paramsab))[1] = 0.5*(sig[0]+sig[1]);   /* sigma   */
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
  L=interval[1]-interval[0];
  FOR_DIM
    h[dim]=PD->h[dim];
  FOR_DIM
    N[dim]=PD->N[dim];

  /* Initialize parallel stuff: fftw + petsc */
  BDD->fft_plan_fw = fftw3d_mpi_create_plan(PETSC_COMM_WORLD,
					    PD->N[2], PD->N[1], PD->N[0],
					    FFTW_FORWARD, FFTW_ESTIMATE);
  BDD->fft_plan_bw = fftw3d_mpi_create_plan(PETSC_COMM_WORLD,
					    PD->N[2], PD->N[1], PD->N[0],
					    FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwnd_mpi_local_sizes(BDD->fft_plan_fw, &local_nx, &local_x_start,
			 &local_ny, &local_y_start, &total_local_size);
  /* Get number of processe */
  MPI_Comm_size(PETSC_COMM_WORLD, &np);

  /* Create Petsc Distributed Array according to fftw data distribution*/
  lz = (PetscInt*) malloc(np*sizeof(*lz));

  MPI_Allgather( &local_nx, 1, MPI_INT, lz, 1, MPI_INT, PETSC_COMM_WORLD);
  ly[0]=PD->N[1];
  lx[0]=PD->N[2];
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     1, 1, np,
	     1,0,
	     lx, ly, lz,
	     &(BDD->da));




  da = BDD->da;
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
  if( verbosity >2)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Subgrids on processes:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "id %d of %d: %d %d %d\t%d %d %d\tfft: %d %d\n",
			      PD->id, PD->np, x[0], x[1], x[2], n[0], n[1], n[2],
			      local_nx, local_x_start);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }
  assert( n[0]*n[1]*n[2] == total_local_size);
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
	    // exp(-beta* Lennard_Jones( r_s, BDD->LJ_paramsa));
	    exp(-beta* Lennard_Jones( r_s, epsilona, sigmaa));
	  gbini_vec[i[2]][i[1]][i[0]] *=
	    // exp(-beta* Lennard_Jones( r_s, BDD->LJ_paramsb));
	    exp(-beta* Lennard_Jones( r_s, epsilonb, sigmab));
	  gabini_vec[i[2]][i[1]][i[0]] *=
	    // exp(-beta* Lennard_Jones( r_s, BDD->LJ_paramsab));
	    exp(-beta* Lennard_Jones( r_s, epsilonab, sigmaab));

	   FOR_DIM
	    {

	      fa_vec[dim][i[2]][i[1]][i[0]] +=
		// Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsa);
		Lennard_Jones_grad( r_s, r[dim], epsilona, sigmaa);
	      fb_vec[dim][i[2]][i[1]][i[0]] +=
		// Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsb);
		Lennard_Jones_grad( r_s, r[dim], epsilonb, sigmab);
	      fab_vec[dim][i[2]][i[1]][i[0]] +=
		// Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsab);
		Lennard_Jones_grad( r_s, r[dim], epsilonab, sigmaab);
	    }

/* 	  FOR_DIM */
/* 	    { */
/* 	      r[dim] = i[dim]*h[dim]; */
/* 	      if( i[dim]>=N[dim]/2) */
/* 		r[dim] -= L; */
/* 	    } */

/* 	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) ); */

/* 	  FOR_DIM */
/* 	    { */

/* 	      fa_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsa); */
/* 	      fb_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsb); */
/* 	      fab_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsab); */
/* 	    } */


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



/*   VecView(BDD->gb_ini,PETSC_VIEWER_STDERR_WORLD);     */
/*   exit(1);   */



  if(BDD->fft_plan_fw == NULL || BDD->fft_plan_bw == NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Failed to get fft_plan of proc %d.\n",
		  PD->id);
      exit(1);
    }


  /* Allocate memory for fft */
  FOR_DIM
    {
      BDD->fg2_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
    }

  BDD->g_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BDD->gfg2_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BDD->fft_scratch = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));





  free(lz);

  return BDD;
}



void BGY3dDiatomicABData_free(BGY3dDiatomicABData BDD)
{
  MPI_Barrier( PETSC_COMM_WORLD);

  FOR_DIM
    {
      VecDestroy(BDD->fa[dim]);
      VecDestroy(BDD->fb[dim]);
      VecDestroy(BDD->fab[dim]);
      VecDestroy(BDD->v[dim]);
      free(BDD->fg2_fft[dim]);
    }
  free(BDD->g_fft);
  free(BDD->gfg2_fft);
  free(BDD->fft_scratch);
  VecDestroy(BDD->ga_ini);
  VecDestroy(BDD->gb_ini);
  VecDestroy(BDD->gab_ini);
  DADestroy(BDD->da);
  // free(BDD->LJ_paramsa);
  // free(BDD->LJ_paramsb);
  // free(BDD->LJ_paramsab);

  fftwnd_mpi_destroy_plan(BDD->fft_plan_fw);
  fftwnd_mpi_destroy_plan(BDD->fft_plan_bw);

  free(BDD);
}


void ComputeDiatomicAB_g_old(BGY3dDiatomicABData BDD, Vec g, Vec g0, Vec dg)
{
  DA da;
  int x[3], n[3], i[3];
  PetscScalar ***g_vec, ***dg_vec;
  // real g_norm;

  da = BDD->da;



  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */

  DAVecGetArray(da, g, &g_vec);
  DAVecGetArray(da, dg, &dg_vec);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  g_vec[i[2]][i[1]][i[0]] = exp(-dg_vec[i[2]][i[1]][i[0]]);
	}
  DAVecRestoreArray(da, g, &g_vec);
  DAVecRestoreArray(da, dg, &dg_vec);

  /* g=g0*exp(-dg) */
  VecPointwiseMult(g, g, g0);

  /* apply same normalization for all g functions */
/*   VecSum(g, &g_norm); */
/*   VecScale(g, 1./(BDD->norm_const*g_norm)); */
  //PetscPrintf(PETSC_COMM_WORLD,"Hier:: %e\n",BDD->norm_const*g_norm);

/*   VecView(g0,PETSC_VIEWER_STDERR_WORLD); */
/*   VecView(g,PETSC_VIEWER_STDERR_WORLD); */
}

void ComputeDiatomicAB_g(Vec g, Vec g0, Vec dg)
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
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  fftw_complex *(fg2_fft[3]), *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, rho, h, sign;

  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BDD->g_fft;
  dg_fft = BDD->gfg2_fft;
  scratch = BDD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  rho = PD->rho;
  fac = L/(2.*M_PI);  /* BDD->f ist nur grad U, nicht F=-grad U  */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* F1*g1a g1b*/
  /************************************************/
  //VecCopy(g1a, dg_help);
  //ShiftVec(da, dg_help, BDD->v[0], N);

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], g1a, f1[dim]);
      if(sign1 != 1.0)
	VecScale(BDD->v[dim], sign1);
      ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, BDD->v[dim], fg2_fft[dim],
			     scratch, x, n, 0);
    }

  /* fft(g-1) */
  //VecCopy(g1b, dg_help);
  //VecShift(dg_help, -1.0);
  ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, g1b, g_fft, scratch,
			 x, n, 0);

/*   PetscPrintf(PETSC_COMM_WORLD, "1: int g= %e\tint fg= %e\n",  */
/* 	      g_fft[0].re*BDD->norm_const,  */
/* 	      fg2_fft[0][0].re*BDD->norm_const); */
  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft[index].re= 0;
	  dg_fft[index].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft[index].re = 0;//g_fft[0].re*h;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = h*h*rho*fac/k;
	      /* phase shift factor for x=x+L/2 */
	      sign = cos(M_PI*ic[0])*cos(M_PI*ic[1])*cos(M_PI*ic[2]);

	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac * sign *
		(fg2_fft[dim][index].re * g_fft[index].im
		 + fg2_fft[dim][index].im * g_fft[index].re) ;

	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac * sign *
		(-fg2_fft[dim][index].re * g_fft[index].re
		 + fg2_fft[dim][index].im * g_fft[index].im);

	    }
	  //fprintf(stderr,"%e\n",dg_fft[index].re);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BDD->fft_plan_bw, dg_help, dg_fft,
			 scratch, x, n, 0.0);

  VecScale(dg_help, PD->beta/L/L/L);

  VecCopy(dg_help, dg);



  /************************************************/
  /* F2*g2a g2b*/
  /************************************************/
  //VecCopy(g2a, dg_help);
  //ShiftVec(da, dg_help, BDD->v[0], N);

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], g2a, f2[dim]);
      if(sign2 != 1.0)
	VecScale(BDD->v[dim], sign2);
      ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, BDD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }

  /* fft(g-1) */
  //VecCopy(g2b, dg_help);
  //VecShift(dg_help, -1.0);
  ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, g2b, g_fft, scratch,
			 x, n, 0);

/*   PetscPrintf(PETSC_COMM_WORLD, "2: int g= %e\tint fg= %e\n",  */
/* 	      g_fft[0].re*BDD->norm_const,  */
/* 	      fg2_fft[0][0].re*BDD->norm_const); */
  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft[index].re= 0;
	  dg_fft[index].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft[index].re = 0;//g_fft[0].re*h;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = h*h*rho*fac/k;
	      /* phase shift factor for x=x+L/2 */
	      sign = cos(M_PI*ic[0])*cos(M_PI*ic[1])*cos(M_PI*ic[2]);

	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac * sign *
		(fg2_fft[dim][index].re * g_fft[index].im
		 + fg2_fft[dim][index].im * g_fft[index].re);

	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac * sign *
		(-fg2_fft[dim][index].re * g_fft[index].re
		 + fg2_fft[dim][index].im * g_fft[index].im);

	    }
	  //fprintf(stderr,"%e\n",dg_fft[index].im);
	  index++;
	}


  ComputeVecfromFFT_fftw(da, BDD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, PD->beta/L/L/L);

/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

  VecAXPY(dg,1.0, dg_help);



  //ShiftVec(da, dg, BDD->v[0], N);
  //VecView(dg,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);

}


/* Compute intramolecular part */
void Compute_dg_Pair_intra(BGY3dDiatomicABData BDD, Vec f[3], Vec g1, Vec g2,
			   Vec dg, Vec dg_help)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  fftw_complex *(fg2_fft[3]), *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, beta; // sign;

  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BDD->g_fft;
  dg_fft = BDD->gfg2_fft;
  scratch = BDD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/
  //VecCopy(g1, dg_help);
  //ShiftVec(da, dg_help, BDD->v[0], N);

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], g1, f[dim]);
      ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, BDD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }

  /* fft(g) */
  //VecCopy(g2, dg_help);
  //VecShift(dg_help, -1.0);
  ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, g2, g_fft, scratch,
			 x, n, 0);

  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft[index].re= 0;
	  dg_fft[index].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft[index].re = 0;//g_fft[0].re*h;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = beta*fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;

	      /* phase shift factor for x=x+L/2 */
	      //sign = cos(M_PI*ic[0])*cos(M_PI*ic[1])*cos(M_PI*ic[2]);



 	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac *
		(h*fg2_fft[dim][index].im * sin(k)/k);

	      /* - should be correct ??? */
	      //dg_fft[index].re -= h*g_fft[index].re*sin(k)/k;



	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac *
		(-h*fg2_fft[dim][index].re * sin(k)/k);



		}
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BDD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, 1./L/L/L);

  //VectorMovetoZero( BDD, dg_help);

  //VecAXPY(dg, 1.0, dg_help);
  VecCopy(dg_help,dg);

  //ShiftVec(da, dg, BDD->v[0], N);
/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}



/* Compute intramolecular part */
void Compute_dg_Pair_intra_ln(BGY3dDiatomicABData BDD, Vec g, real sign, Vec dg, Vec dg_help)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3], local_size;
  fftw_complex *g_fft, *dg_fft, *scratch;
  real fac, L, k, h, beta;
  PetscScalar *g_vec;


  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BDD->g_fft;
  dg_fft = BDD->gfg2_fft;
  scratch = BDD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  fac = L/(2.*M_PI); /* siehe oben ... */
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* F(g) */
  /************************************************/


  ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, g, g_fft, scratch,
			 x, n, 0);

  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft[index].re= 0;
	  dg_fft[index].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft[index].re = g_fft[0].re*h;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      //k_fac = beta*fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;

	      /* + should be correct ??? */
	      dg_fft[index].re += h*g_fft[index].re*sin(k)/k;

		}
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BDD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, 1./L/L/L);

  /* ln(g) */
  VecGetArray( dg_help, &g_vec);
  VecGetLocalSize(dg_help, &local_size);

  for(index=0; index<local_size; index++)
    g_vec[index] = sign*log(g_vec[index]);
  //PetscPrintf(PETSC_COMM_WORLD, "%e\n", g_vec[0]);
  VecRestoreArray(dg_help, &g_vec);
  /******************************/

  VecAXPY(dg, 1.0, dg_help);
  //VecCopy(dg_help,dg);


/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}




/* Compute normalization condition */
void Compute_dg_Pair_normalization_intra(BGY3dDiatomicABData BDD, Vec g,
					 Vec dg, Vec dg_help)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  fftw_complex  *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, beta; // sign;

  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];


  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BDD->g_fft;
  dg_fft = BDD->gfg2_fft;
  scratch = BDD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* fft(g/t) */
  ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, g, g_fft, scratch,
			 x, n, 0);



  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft[index].re= 0;
	  dg_fft[index].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {

	      dg_fft[index].re = h*g_fft[0].re;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = beta*fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;



	      /* + oder - hier ??? */
	      dg_fft[index].re += h*g_fft[index].re*sin(k)/k;



	    }
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BDD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);


  VecScale(dg_help, 1./L/L/L);


  //VecAXPY(dg, 1.0, dg_help);
  VecCopy(dg_help,dg);


  //VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);

}

/* Compute normalization condition */
void Compute_dg_Pair_normalization(BGY3dDiatomicABData BDD, Vec g1, Vec g2,
				   Vec dg, Vec dg_help)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  fftw_complex  *(fg2_fft[3]), *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, sign, beta;

  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BDD->g_fft;
  dg_fft = BDD->gfg2_fft;
  scratch = BDD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* fft(g) */
  //VecCopy(g2, dg_help);
  //VecShift(dg_help, -1.0);
  ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, g1, g_fft, scratch,
			 x, n, 0);

  ComputeFFTfromVec_fftw(da, BDD->fft_plan_fw, g2, fg2_fft[0], scratch,
			 x, n, 0);

  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft[index].re= 0;
	  dg_fft[index].im= 0;

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft[index].re = h*h*g_fft[0].re*fg2_fft[0][0].re;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = beta*fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;
	      /* phase shift factor for x=x+L/2 */
	      sign = cos(M_PI*ic[0])*cos(M_PI*ic[1])*cos(M_PI*ic[2]);

	      dg_fft[index].re += h*h*sign*
		(g_fft[index].re*fg2_fft[0][index].re -
		 g_fft[index].im*fg2_fft[0][index].im );

	      dg_fft[index].im += h*h*sign*
		(g_fft[index].im*fg2_fft[0][index].re +
		 g_fft[index].re*fg2_fft[0][index].im );

	    }
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BDD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, 1./L/L/L/L/L/L);


  //VecAXPY(dg, 1.0, dg_help);
  VecCopy(dg_help,dg);


  //VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);

}


void Solve_Normalization_old(BGY3dDiatomicABData BDD, Vec ga, Vec gb, Vec gab, Vec gba,
			 Vec ta, Vec tb, Vec tab, Vec tba, Vec dg, Vec dg_help,
			 real norm_tol)
{
  int i;
  real tab_norm, ta_norm, tb_norm, h; // a=0.2, g_norm, h;

  h = BDD->PD->h[0]*BDD->PD->h[1]*BDD->PD->h[2]/1000.;


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


void Solve_Normalization(BGY3dDiatomicABData BDD, Vec gc, Vec g, Vec t, Vec dg,
			     Vec dg_help)
{

  //VecCopy(g, t);


  Compute_dg_Pair_normalization_intra( BDD, gc, dg, dg_help);
  VecPointwiseDivide(t, g, dg);

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
Vec BGY3d_solve_DiatomicAB(ProblemData *PD, Vec g_ini, int vdim)
{
  BGY3dDiatomicABData BDD;
  Vec g0a, g0b, g0ab, dga, dgb, dgab, dg_new, dg_new2, f, ga, gb, gab;
  Vec dgba, gba;
  Vec ta, tb, tab, tba;
  real a=0.9, h;
  int max_iter=25, iter;
  PetscScalar dga_norm, dgb_norm, dgab_norm, norm_tol=1.0e-6;

  PetscScalar dgba_norm;
  PetscTruth kflg; //, load_flag;
  PetscViewer viewer;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving Molecular BGY3d equation with Fourier ansatz...\n");

  kflg = bgy3d_getopt_test ("-pair");
  if(kflg)
    {

      BDD = BGY3dDiatomicABData_Pair_malloc(PD);
    }
  else
    {
      exit(1);
      //BDD = BGY3dFourierData_malloc(PD);
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }

  /* read BGY3dDiv specific things from command line */
  /* Mixing parameter */
  bgy3d_getopt_real ("-lambda", &a);
  if(a>1 || a<0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lambda out of range: lambda=%f\n",a);
      exit(1);
    }

  /* Number of total iterations */
  bgy3d_getopt_int ("-max_iter", &max_iter);
  /* norm_tol for convergence test */
  bgy3d_getopt_real ("-norm_tol", &norm_tol);

  /*********************************/

  PetscPrintf(PETSC_COMM_WORLD,"lambda = %f\n",a);
  PetscPrintf(PETSC_COMM_WORLD,"tolerance = %e\n",norm_tol);
  PetscPrintf(PETSC_COMM_WORLD,"max_iter = %d\n",max_iter);

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

/*   PetscOptionsHasName(PETSC_NULL,"-loadpair",&load_flag); */
/*   if(load_flag) */
/*     { */
/*       PetscViewerBinaryOpen(PETSC_COMM_WORLD,"veca.m", */
/* 			    FILE_MODE_READ,&viewer); */
/*       VecLoad(viewer,VECMPI, &dgb); */
/*       PetscViewerDestroy(viewer); */
/*       PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecb.m", */
/* 			    FILE_MODE_READ,&viewer); */
/*       VecLoad(viewer,VECMPI, &dga); */
/*       PetscViewerDestroy(viewer); */
/*       PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m", */
/* 			    FILE_MODE_READ,&viewer); */
/*       VecLoad(viewer,VECMPI, &dgab); */
/*       PetscViewerDestroy(viewer); */
/*     } */


/*   VecView(v,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */


/*   VecSum(g0a, &dga_norm); */
/*   VecScale(g0a, 1./(BDD->norm_const*dga_norm)); */
/*   VecSum(g0b, &dgb_norm); */
/*   VecScale(g0b, 1./(BDD->norm_const*dgb_norm)); */
/*   VecSum(g0ab, &dgab_norm); */
/*   VecScale(g0ab, 1./(BDD->norm_const*dgab_norm)); */
//  BDD->c_ab = 0.95;//dga_norm/dgb_norm;
//  BDD->c_aab = 0.98;//dga_norm/dgab_norm;

  h = PD->h[0]*PD->h[1]*PD->h[2]/1000.;

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
      if(kflg)
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

      if(dga_norm<=norm_tol &&  dgb_norm<=norm_tol && dgab_norm<=norm_tol &&
	 dgba_norm<=norm_tol)
	break;

/*       VecSum(ga, &dga_norm); */
/*       VecSum(gb, &dgb_norm); */
/*       VecSum(gab, &dgab_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD," %e %e %e ll %e %e\n",  */
/* 		  dga_norm*h, dgb_norm*h, dgab_norm*h, */
/* 		  dga_norm/dgb_norm, dgb_norm/dga_norm); */



    }
/*   VecSet(dg_new,0.0); */
/*   Compute_dg_Pair_intra(BDD, BDD->fa, ga, gb, dg_new, f); */
/*   VecCopy(dg_new, gab); */
/*   VecSet(dg_new,0.0); */
/*   Compute_dg_Pair_intra(BDD, BDD->fb, gb, ga, dg_new, f); */
/*   VecCopy(dg_new, gba); */


/*   Compute_dg_Pair_inter(BDD, BDD->fa, ga, gba, BDD->fab, gba, gb,  */
/* 	  			dgab, f); */
/*   ComputeDiatomicAB_g(BDD, gab, g0ab, dgab); */
/*   Compute_dg_Pair_inter(BDD, BDD->fb, gb, gab, BDD->fab, gab, ga,  */
/* 				dgba, f); */
/*   ComputeDiatomicAB_g(BDD, gba, g0ab, dgba); */
/*   Compute_dg_Pair_intra(BDD, BDD->fab, gab, gab, dg_new, f); */
/*   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecdg.m",&viewer); */
/*   PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); */
/*   VecView(dg_new,viewer); */
/*   PetscViewerDestroy(viewer); */

  /* g=g0*exp(-dg) */
/*   ComputeDiatomicAB_g(BDD, ga, g0a, dga); */
/*   ComputeDiatomicAB_g(BDD, gb, g0b, dgb); */
/*   ComputeDiatomicAB_g(BDD, gab, g0ab, dgab); */
  //VecCopy(dg, g);
/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
  if( 0 && kflg)
    {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2_new.bin",
			    FILE_MODE_WRITE,&viewer);
      VecCopy(ga,dga);
      //ShiftVec(BDD->da, dga, BDD->v[0], PD->N);
      VecView(dga,viewer);
      PetscViewerDestroy(viewer);
      PetscPrintf(PETSC_COMM_WORLD,"g2 vector written to file \"g2_new.bin\".\n");
    }
  /*************************************/
  /* output */
  /* g_a */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"veca.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"veca.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(ga,viewer);
  PetscViewerDestroy(viewer);
  /* g_b */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecb.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecb.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gb,viewer);
  PetscViewerDestroy(viewer);
  /* g_ab */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecab.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gab,viewer);
  PetscViewerDestroy(viewer);
  /* g_ba */
  //PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecba.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m",FILE_MODE_WRITE,&viewer);
  //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  //VecView(gba,viewer);
  //PetscViewerDestroy(viewer);
  /************************************/
/*   Compute_dg_Pair_normalization(BDD, ga, gab, dg_new, f); */
/*   PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vec.m",&viewer); */
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m",FILE_MODE_WRITE,&viewer);
/*   PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); */
/*   VecView(dg_new,viewer); */
/*   PetscViewerDestroy(viewer); */


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

fftw_complex *ComputeFFTfromVec_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
				fftw_complex *g_fft, fftw_complex *scratch,
				int x[3], int n[3], real c)
{
  int index=0, i[3];
  PetscScalar ***g_vec;

  if(g_fft==NULL)
    g_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(*g_fft));

  DAVecGetArray(da, g, &g_vec);
  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index=0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  //g_fft[index].re = g_vec[i[2]][i[1]][i[0]]/2.+c;
	  g_fft[index].re = g_vec[i[2]][i[1]][i[0]]+c;
	  g_fft[index].im = 0;
	  index++;
	}
  DAVecRestoreArray(da, g, &g_vec);
  /* forward fft */
  fftwnd_mpi( fft_plan, 1, g_fft, scratch, FFTW_NORMAL_ORDER);

  return g_fft;

}


void ComputeVecfromFFT_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
			    fftw_complex *g_fft, fftw_complex *scratch,
			    int x[3], int n[3], real c)
{
  int index=0, i[3];
  PetscScalar ***g_vec;

  if(g_fft==NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Error: g_fft==NULL!\n");
      exit(1);
    }

  /* backward fft */
  fftwnd_mpi( fft_plan, 1, g_fft, scratch, FFTW_NORMAL_ORDER);

  DAVecGetArray(da, g, &g_vec);
  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index=0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  //g_vec[i[2]][i[1]][i[0]] = g_fft[index].re*2.+c; // Factor 2!!!
	  g_vec[i[2]][i[1]][i[0]] = g_fft[index].re+c;
	  index++;
	}
  DAVecRestoreArray(da, g, &g_vec);


}


