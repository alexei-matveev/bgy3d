/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.6 2007-01-12 10:20:43 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"


#define rab 1.0



BGY3dDiatomicABData BGY3dDiatomicABData_Pair_malloc(PData PD)
{
  BGY3dDiatomicABData BDD;
  DA da;
  real interval[2], h[3], N[3], L, r[3], r_s, beta;
  int i[3], x[3], n[3], dim;
  PetscScalar ***(fa_vec[3]),***(fb_vec[3]),***(fab_vec[3]);
  PetscScalar ***gaini_vec, ***gbini_vec, ***gabini_vec;
  int bufsize, np;
  int local_nx, local_x_start, local_ny, local_y_start, total_local_size;
  PetscInt *lx, ly[1], lz[1];


  real eps[]={ 1.0, 1.0}, sig[]={ 0.5, 1.0}; 
  //real eps[]={ 1.0, 1.0}, sig[]={ 1.0, 1.0}; 


  BDD = (BGY3dDiatomicABData) malloc(sizeof(*BDD));
 

  BDD->LJ_paramsa = (void* ) malloc(sizeof(real)*2);
  ((real*)(BDD->LJ_paramsa))[0] = eps[0];   /* espilon */
  ((real*)(BDD->LJ_paramsa))[1] = sig[0];   /* sigma   */

  BDD->LJ_paramsb = (void* ) malloc(sizeof(real)*2);
  ((real*)(BDD->LJ_paramsb))[0] = eps[1];   /* espilon */
  ((real*)(BDD->LJ_paramsb))[1] = sig[1];   /* sigma   */

  BDD->LJ_paramsab = (void* ) malloc(sizeof(real)*2);
  ((real*)(BDD->LJ_paramsab))[0] = sqrt(eps[0]*eps[1]);   /* espilon */
  ((real*)(BDD->LJ_paramsab))[1] = 0.5*(sig[0]+sig[1]);   /* sigma   */



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
					    PD->N[0], PD->N[1], PD->N[2],
					    FFTW_FORWARD, FFTW_ESTIMATE);
  BDD->fft_plan_bw = fftw3d_mpi_create_plan(PETSC_COMM_WORLD, 
					    PD->N[1], PD->N[0], PD->N[2],
					    FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwnd_mpi_local_sizes(BDD->fft_plan_fw, &local_nx, &local_x_start,
			 &local_ny, &local_y_start, &total_local_size);
  /* Get number of processe */
  MPI_Comm_size(PETSC_COMM_WORLD, &np);
  
  /* Create Petsc Distributet Array according to fftw data distribution*/
  lx = (PetscInt*) malloc(np*sizeof(*lx));
  MPI_Allgather( &local_nx, 1, MPI_INT, lx, np, PETSC_COMM_WORLD);
  ly[0]=PD->N[1];
  lz[0]=PD->N[2];
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2], 
	     np, 1, 1,
	     1,1,
	     lx, ly, lz,
	     &(BDD->da));
  free(lx);



  da = BDD->da;
  
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
  

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
  
  if( verbosity >2)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Subgrids on processes:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "id %d of %d: %d %d %d\t%d %d %d\n", 
			      PD->id, PD->np, x[0], x[1], x[2], n[0], n[1], n[2]);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
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
	    exp(-beta* Lennard_Jones( r_s, BDD->LJ_paramsa));
	  gbini_vec[i[2]][i[1]][i[0]] *= 
	    exp(-beta* Lennard_Jones( r_s, BDD->LJ_paramsb));
	  gabini_vec[i[2]][i[1]][i[0]] *= 
	    exp(-beta* Lennard_Jones( r_s, BDD->LJ_paramsab));

	   FOR_DIM
	    {
	      
	      fa_vec[dim][i[2]][i[1]][i[0]] += 
		Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsa);
	      fb_vec[dim][i[2]][i[1]][i[0]] += 
		Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsb);
	      fab_vec[dim][i[2]][i[1]][i[0]] += 
		Lennard_Jones_grad( r_s, r[dim], BDD->LJ_paramsab);
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
  

  /* Create plan for 3d fft */
  BDD->fft_plan = fft_3d_create_plan(PETSC_COMM_WORLD,
				     PD->N[0], PD->N[1], PD->N[2],
				     x[0], x[0]+n[0]-1,
				     x[1], x[1]+n[1]-1,
				     x[2], x[2]+n[2]-1,
				     x[0], x[0]+n[0]-1,
				     x[1], x[1]+n[1]-1,
				     x[2], x[2]+n[2]-1,
				     0,
				     0,
				     &bufsize);
  if(BDD->fft_plan == NULL) 
    {
      PetscPrintf(PETSC_COMM_WORLD, "Failed to get fft_plan of proc %d.\n",
		  PD->id);
      exit(1);
    }
  
  
  /* Allocate memory for fft */
  FOR_DIM
    {
      BDD->fg2_fft[dim] = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));     
    }
  
  BDD->g_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));
  BDD->gfg2_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));

  


  


  return BDD;
}
  	  


void BGY3dDiatomicABData_free(BGY3dDiatomicABData BDD)
{
  int dim;
  
 
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
  VecDestroy(BDD->ga_ini);
  VecDestroy(BDD->gb_ini);
  VecDestroy(BDD->gab_ini);
  DADestroy(BDD->da);
  free(BDD->LJ_paramsa);
  free(BDD->LJ_paramsb);
  free(BDD->LJ_paramsab);

  fft_3d_destroy_plan(BDD->fft_plan); 

  free(BDD);
}
  

void ComputeDiatomicAB_g(BGY3dDiatomicABData BDD, Vec g, Vec g0, Vec dg)
{
  DA da;
  int x[3], n[3], i[3];
  PetscScalar ***g_vec, ***dg_vec;
  real g_norm;

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

void Compute_dg_Pair_inter(BGY3dDiatomicABData BDD, Vec f1[3], Vec g1a, Vec g1b, 
			   Vec f2[3], Vec g2a, Vec g2b, Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3];
  FFT_DATA *(fg2_fft[3]), *g_fft, *dg_fft;
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
  L = PD->interval[1]-PD->interval[0];
  rho = PD->rho;
  fac = -L/(2.*M_PI);


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
      ComputeFFTfromVec(da, BDD->fft_plan, BDD->v[dim], fg2_fft[dim], x, n, 0);
    }
  
  /* fft(g-1) */
  //VecCopy(g1b, dg_help); 
  //VecShift(dg_help, -1.0);
  ComputeFFTfromVec(da, BDD->fft_plan, g1b, g_fft, x, n, 0);
  
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
  ComputeVecfromFFT(da, BDD->fft_plan, dg_help, dg_fft, 
		    x, n, 0.0);
  
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
      ComputeFFTfromVec(da, BDD->fft_plan, BDD->v[dim], fg2_fft[dim], x, n, 0);
    }
  
  /* fft(g-1) */
  //VecCopy(g2b, dg_help); 
  //VecShift(dg_help, -1.0);
  ComputeFFTfromVec(da, BDD->fft_plan, g2b, g_fft, x, n, 0);
  
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
	  //fprintf(stderr,"%e\n",dg_fft[index].re);
	  index++;
	}


  ComputeVecfromFFT(da, BDD->fft_plan, dg_help, dg_fft, 
		    x, n, 0.0);

  VecScale(dg_help, PD->beta/L/L/L);

/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */
 
  VecAXPY(dg,1.0, dg_help);
  //VecWAXPY(dg, 1.0, dg,dg_help);


  //ShiftVec(da, dg, BDD->v[0], N);
/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}
  

/* Compute intramolecular part */
void Compute_dg_Pair_intra(BGY3dDiatomicABData BDD, Vec f[3], Vec g1, Vec g2, 
			   Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3];
  FFT_DATA *(fg2_fft[3]), *g_fft, *dg_fft;
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
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = -L/(2.*M_PI);
  

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
      ComputeFFTfromVec(da, BDD->fft_plan, BDD->v[dim], fg2_fft[dim], x, n, 0);
    }
  
  /* fft(g) */
  //VecCopy(g2, dg_help); 
 
  ComputeFFTfromVec(da, BDD->fft_plan, g2, g_fft, x, n, 0);
  
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
	     
	      /* + oder - hier ??? */
	      dg_fft[index].re += h*g_fft[index].re*sin(k)/k;
	      
	      

	      FOR_DIM 
		dg_fft[index].im += ic[dim] * k_fac *
		(-h*fg2_fft[dim][index].re * sin(k)/k);
		 

	      
		}
	  //fprintf(stderr,"%e\n",dg_fft[index].re);
	  index++;
	}
  ComputeVecfromFFT(da, BDD->fft_plan, dg_help, dg_fft, 
		    x, n, 0.0);
  
  VecScale(dg_help, 1./L/L/L);

  VecAXPY(dg, 1.0, dg_help); 


  //ShiftVec(da, dg, BDD->v[0], N);
/*   VecView(BDD->v[0],PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

/* solve with product ansatz g=g0*dg */
Vec BGY3d_solve_DiatomicAB(PData PD, Vec g_ini, int vdim)
{
  BGY3dDiatomicABData BDD;
  Vec g0a, g0b, g0ab, dga, dgb, dgab, dg_new, dg_new2, f, ga, gb, gab;
  real a=0.9, h;
  int max_iter=25, iter;
  PetscScalar dga_norm, dgb_norm, dgab_norm, norm_tol=1.0e-6;
  PetscTruth kflg, load_flag;
  PetscViewer viewer;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving Molecular BGY3d equation with Fourier ansatz...\n");
  
  PetscOptionsHasName(PETSC_NULL,"-pair",&kflg);
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
  PetscOptionsGetReal(PETSC_NULL,"-lambda",&a, PETSC_NULL);
  if(a>1 || a<0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lambda out of range: lambda=%f\n",a);
      exit(1);
    }
  
  /* Number of total iterations */
  PetscOptionsGetInt(PETSC_NULL,"-max_iter",&max_iter, PETSC_NULL);
  /* norm_tol for convergence test */
  PetscOptionsGetReal(PETSC_NULL,"-norm_tol",&norm_tol, PETSC_NULL);

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
  

//  VecSum(g0a, &dga_norm);
//  VecSum(g0b, &dgb_norm);
//   VecSum(g0ab, &dgab_norm);
//  BDD->c_ab = 0.95;//dga_norm/dgb_norm;
//  BDD->c_aab = 0.98;//dga_norm/dgab_norm;

  h = PD->h[0]*PD->h[1]*PD->h[2]/1000.;

  /* g=g0+exp(-dg) */
  ComputeDiatomicAB_g(BDD, ga, g0a, dga);
  ComputeDiatomicAB_g(BDD, gb, g0b, dgb);
  ComputeDiatomicAB_g(BDD, gab, g0ab, dgab);
  
/*   VecScale(gb, BDD->c_ab); */
/*   VecScale(gab, BDD->c_aab); */

  for(iter=0; iter<max_iter; iter++)
    {
     
      
     
      /* f=integral(g) */
      if(kflg)
	{

	   /* g_ab */
	  Compute_dg_Pair_inter(BDD, BDD->fa, ga, gab, BDD->fab, gab, gb, 
				dg_new, f);
	  Compute_dg_Pair_intra(BDD, BDD->fa, ga, gb, dg_new, f);
	  Compute_dg_Pair_inter(BDD, BDD->fb, gb, gab, BDD->fab, gab, ga,  
 				dg_new2, f);
	  Compute_dg_Pair_intra(BDD, BDD->fb, gb, ga, dg_new2, f);
	  VecAXPBY(dg_new, 0.5, 0.5, dg_new2);
	  VecWAXPY(f, -1.0, dg_new, dgab);
	  VecNorm(f, NORM_INFINITY, &dgab_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms: ab= %e  ",  
		      iter+1, dgab_norm);
	  VecAXPBY(dgab, a, (1-a), dg_new);
	  ComputeDiatomicAB_g(BDD, gab, g0ab, dgab);
	  //VecScale(gab, BDD->c_aab);
	   /* g_ba */
/* 	  Compute_dg_Pair_inter(BDD, BDD->fb, gb, gab, BDD->fab, gab, ga,  */
/* 				dg_new, f); */
/* 	  Compute_dg_Pair_intra(BDD, BDD->fb, gb, ga, dg_new, f); */
/* 	  VecWAXPY(f, -1.0, dg_new, dgab); */
/* 	  VecNorm(f, NORM_INFINITY, &dgab_norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms: ab= %e  ",   */
/* 		      iter+1, dgab_norm); */
/* 	  VecAXPBY(dgab, a, (1-a), dg_new); */
/* 	  ComputeDiatomicAB_g(BDD, gab, g0ab, dgab); */
	  /* g_a */
	  Compute_dg_Pair_inter(BDD, BDD->fa, ga, ga, BDD->fab, gab, gab, 
				dg_new, f);
	  Compute_dg_Pair_intra(BDD, BDD->fab, gab, gab, dg_new, f);
	  VecWAXPY(f, -1.0, dg_new, dga);
	  VecNorm(f, NORM_INFINITY, &dga_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"a= %e  ", dga_norm);
	  VecAXPBY(dga, a, (1-a), dg_new);
	  ComputeDiatomicAB_g(BDD, ga, g0a, dga);
	  /* g_b */
	  Compute_dg_Pair_inter(BDD, BDD->fb, gb, gb, BDD->fab, gab, gab, 
				dg_new, f);
	   Compute_dg_Pair_intra(BDD, BDD->fab, gab, gab, dg_new, f);
	  VecWAXPY(f, -1.0, dg_new, dgb);
	  VecNorm(f, NORM_INFINITY, &dgb_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"b= %e\n", dgb_norm);
	  VecAXPBY(dgb, a, (1-a), dg_new);
	  ComputeDiatomicAB_g(BDD, gb, g0b, dgb);
	  //VecScale(gb, BDD->c_ab);
	}
      else
	;
	
      
      

      if(dga_norm<=norm_tol && dgb_norm<=norm_tol && dgab_norm<=norm_tol)
	break;
   
/*       VecSum(ga, &dga_norm); */
/*       VecSum(gb, &dgb_norm); */
/*       VecSum(gab, &dgab_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD," %e %e %e ll %e %e\n",  */
/* 		  dga_norm*h, dgb_norm*h, dgab_norm*h, */
/* 		  dga_norm/dgb_norm, dgb_norm/dga_norm); */
     
      

    }

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
  /************************************/


  VecDestroy(ga);
  VecDestroy(gb);
  VecDestroy(gab);
  VecDestroy(dga);
  VecDestroy(dgb);
  VecDestroy(dgab);
  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(f);
  
  // ExtractAxis(BDD, g, 0);


  BGY3dDiatomicABData_free(BDD);
  

  return PETSC_NULL;
}
