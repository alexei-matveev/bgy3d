/*==========================================================*/
/*  $Id: bgy3dfourier.c,v 1.14 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fft.h"
#include "bgy3ddiv.h"
#include "bgy3dfourier.h"

BGY3dFourierData BGY3dFourierData_malloc(ProblemData *PD)
{
  BGY3dFourierData BDD;
  DA da;
  real interval[2], h[3], N[3], L, r[3], r_s, **x_M, beta; //, h_g2;
  int i[3], x[3], n[3], N_M, k,  N_g2, index;
  PetscScalar  ***(f_vec[3]), ***gini_vec, ***(v2_vec[3]);
  PetscScalar *g2_vec;
  Vec g2, g2_3d;
  PetscViewer pview;
  int bufsize;
  PetscTruth g2_flg;
  real epsilon, sigma;

  BDD = (BGY3dFourierData) malloc(sizeof(*BDD));

  BDD->LJ_params[0] = 1.0;
  BDD->LJ_params[1] = 1.0;
  epsilon = BDD->LJ_params[0];
  sigma = BDD->LJ_params[1];

  BDD->beta = PD->beta;
  BDD->rho  = PD->rho;
  beta = PD->beta;

  BDD->PD = PD;

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];
  L=interval[1]-interval[0];
  FOR_DIM
    h[dim]=PD->h[dim];
  FOR_DIM
    N[dim]=PD->N[dim];

  /* Create distributed array */
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	     1,1, PETSC_NULL,PETSC_NULL,PETSC_NULL, &(BDD->da));

  da = BDD->da;

  /* Create global vectors */
  DACreateGlobalVector(da, &(BDD->g_ini));

  FOR_DIM
    {
      VecDuplicate(BDD->g_ini, &(BDD->f[dim]));
      VecDuplicate(BDD->g_ini, &(BDD->v[dim]));
    }


  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if( verbosity >2)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Subgrids on processes:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "id %d of %d: %d %d %d\t%d %d %d\n",
			      PD->id, PD->np, x[0], x[1], x[2], n[0], n[1], n[2]);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }


  /* Load molecule from file */
  x_M = Load_Molecule(&N_M);

  g2_flg = bgy3d_getopt_test ("-g2load3d");
  if(g2_flg)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Loading 3d g2 vector.\n");
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2.bin",
			    FILE_MODE_READ,&pview);
      VecLoad(pview,VECMPI, &g2_3d);
      PetscViewerDestroy(pview);

      VecGetSize(g2_3d, &N_g2);
      if( N_g2!=PD->N3)
	{
	  PetscPrintf(PETSC_COMM_WORLD,"Wrong size for g2 vector!\nAborting.\n");
	  exit(1);
	}
    }
  else
    {
      /* Load g2 from file, sequential only */
      /* g2 has to be on a grid [0,L] with L=interval[1]-interval[0] */
      PetscViewerBinaryOpen(PETSC_COMM_SELF,"g2file.bin" , FILE_MODE_READ,
			    &pview);
      VecLoad( pview, VECSEQ, &g2);
      PetscViewerDestroy(pview);
      VecGetSize(g2,&N_g2);
      // h_g2 = L/N_g2;

      //g2_vec[0] = g2_vec[1]+g2_vec[1]-g2_vec[2];
    }



  FOR_DIM
    VecSet(BDD->f[dim],0.0);


  VecSet(BDD->g_ini, 1.0);



  DAVecGetArray(da, BDD->g_ini, &gini_vec);

  FOR_DIM
    {
      DAVecGetArray(da, BDD->f[dim], &(f_vec[dim]));
      DAVecGetArray(da, BDD->v[dim], &(v2_vec[dim]));
    }

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */
	  /* loop over particles and grid */
	  for(k=0; k<N_M; k++)
	    {
	      /* 	      for(index=0; index<9; index++) */
	      /* 		{ */
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim]+interval[0]-x_M[k][dim];
		}



	      /* 	      switch(index) */
	      /* 		{ */
	      /* 		case 0: r[0]+=10; */
	      /* 		  break; */
	      /* 		case 1:  */
	      /* 		  break; */
	      /* 		case 2: r[0]-=10; */
	      /* 		  break; */
	      /* 		case 3: r[0]+=10; r[1]+=10; */
	      /* 		  break; */
	      /* 		case 4:           r[1]+=10;  */
	      /* 		  break;  */
	      /* 		case 5: r[0]-=10; r[1]+=10; */
	      /* 		  break; */
	      /* 		case 6: r[0]+=10; r[1]-=10; */
	      /* 		  break; */
	      /* 		case 7:           r[1]-=10;  */
	      /* 		  break;  */
	      /* 		case 8: r[0]-=10; r[1]-=10; */
	      /* 		  break; */
	      /* 		} */
	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	      FOR_DIM
		f_vec[dim][i[2]][i[1]][i[0]] +=
		Lennard_Jones_grad( r_s, r[dim], epsilon, sigma);




	      gini_vec[i[2]][i[1]][i[0]] *=
		exp(-beta* Lennard_Jones( r_s, epsilon, sigma));
	    }

	  /*	}*/


	}
  if( g2_flg)
    {

      /* loop over local portion of grid */
      for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
	for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
	  for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	    {

	      /* f*g2 */
	      /* this vector lives on a different grid -> [0,L]^3 for fft */
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim];
		  if( i[dim]>=N[dim]/2)
		    r[dim] -= L;
		}

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	      FOR_DIM
		v2_vec[dim][i[2]][i[1]][i[0]] =
		Lennard_Jones_grad( r_s, r[dim], epsilon, sigma);

	    }

      FOR_DIM
	DAVecRestoreArray(da, BDD->v[dim], &(v2_vec[dim]));

      FOR_DIM
	VecPointwiseMult(BDD->v[dim], BDD->v[dim], g2_3d);

      //VecView(g2_3d,PETSC_VIEWER_STDERR_WORLD);
    }
  else
    {
      // only used in this branch, declare inside the block:
      real h_g2 = L / N_g2;

      VecGetArray(g2, &g2_vec);
      /* loop over local portion of grid */
      for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
	for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
	  for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	    {

	      /* f*g2 */
	      /* this vector lives on a different grid -> [0,L]^3 for fft */
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim];
		  if( i[dim]>=N[dim]/2)
		    r[dim] -= L;
		}

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	      index =(int) floor(r_s/h_g2);
	      /* g2=1 if r_s>l */
	      if(index >= N_g2-1)
	    {
	      FOR_DIM
		v2_vec[dim][i[2]][i[1]][i[0]] =
		Lennard_Jones_grad( r_s, r[dim], epsilon, sigma);
	    }
	      else
		{
		  FOR_DIM
		    v2_vec[dim][i[2]][i[1]][i[0]] =
		    Lennard_Jones_grad( r_s, r[dim], epsilon, sigma)*
		    (g2_vec[index]+ (g2_vec[index+1]-g2_vec[index])/h_g2*
		     (r_s-index*h_g2));
		}

	    }
      VecRestoreArray(g2, &g2_vec);
      FOR_DIM
	DAVecRestoreArray(da, BDD->v[dim], &(v2_vec[dim]));
    }



  DAVecRestoreArray(da, BDD->g_ini, &gini_vec);
  FOR_DIM
    {
      DAVecRestoreArray(da, BDD->f[dim], &(f_vec[dim]));
    }




/*   VecView(g2,PETSC_VIEWER_STDERR_WORLD);      */
/*   exit(1);    */




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


  /* Compute fft for f*g2 */
  FOR_DIM
    {
      BDD->fg2_fft[dim] = NULL;
      BDD->fg2_fft[dim] = ComputeFFTfromVec(da, BDD->fft_plan,
					    BDD->v[dim],
					    BDD->fg2_fft[dim]);
    }
  BDD->g_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));
  BDD->gfg2_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));



  Molecule_free(x_M, N_M);
  if(g2_flg)
    VecDestroy(g2_3d);
  else
    VecDestroy(g2);

  return BDD;
}


BGY3dFourierData BGY3dFourierData_kirk_malloc(ProblemData *PD)
{
  BGY3dFourierData BDD;
  DA da;
  real interval[2], h[3], N[3], L, r[3], r_s, beta;
  int i[3], x[3], n[3];
  PetscScalar ***(f_vec[3]), ***gini_vec;
  int bufsize;
  real epsilon, sigma;

  BDD = (BGY3dFourierData) malloc(sizeof(*BDD));

  BDD->LJ_params[0] = 1.0;
  BDD->LJ_params[1] = 1.0;
  epsilon = BDD->LJ_params[0];
  sigma = BDD->LJ_params[1];


  BDD->beta = PD->beta;
  BDD->rho  = PD->rho;
  beta = PD->beta;

  BDD->PD = PD;

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];
  L=interval[1]-interval[0];
  FOR_DIM
    h[dim]=PD->h[dim];
  FOR_DIM
    N[dim]=PD->N[dim];

  /* Create distributed array */
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	     1,1, PETSC_NULL,PETSC_NULL,PETSC_NULL, &(BDD->da));

  da = BDD->da;

  /* Create global vectors */
  DACreateGlobalVector(da, &(BDD->g_ini));
  FOR_DIM
    {
      VecDuplicate(BDD->g_ini, &(BDD->f[dim]));
      VecDuplicate(BDD->g_ini, &(BDD->v[dim]));
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
    VecSet(BDD->f[dim],0.0);


  VecSet(BDD->g_ini, 1.0);



  DAVecGetArray(da, BDD->g_ini, &gini_vec);
  FOR_DIM
    {
      DAVecGetArray(da, BDD->f[dim], &(f_vec[dim]));
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
	  gini_vec[i[2]][i[1]][i[0]] *=
	    exp(-beta* Lennard_Jones( r_s, epsilon, sigma));



	  FOR_DIM
	    {
	      r[dim] = i[dim]*h[dim];
	      if( i[dim]>=N[dim]/2)
		r[dim] -= L;
	    }

	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  FOR_DIM
	    f_vec[dim][i[2]][i[1]][i[0]] +=
	    Lennard_Jones_grad( r_s, r[dim], epsilon, sigma);




	}




  DAVecRestoreArray(da, BDD->g_ini, &gini_vec);

  FOR_DIM
    {
      DAVecRestoreArray(da, BDD->f[dim], &(f_vec[dim]));
    }



/*   VecView(BDD->f[0],PETSC_VIEWER_STDERR_WORLD);  */
/*   VecView(BDD->g_ini,PETSC_VIEWER_STDERR_WORLD);    */
/*   exit(1);  */


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



void BGY3dFourierData_free(BGY3dFourierData BDD)
{
  FOR_DIM
    {
      VecDestroy(BDD->f[dim]);
      VecDestroy(BDD->v[dim]);
      free(BDD->fg2_fft[dim]);
    }
  free(BDD->g_fft);
  free(BDD->gfg2_fft);
  VecDestroy(BDD->g_ini);
  DADestroy(BDD->da);

  fft_3d_destroy_plan(BDD->fft_plan);

  free(BDD);
}


void ExtractAxis(BGY3dFourierData BDD, Vec g, int axis)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], N[3], ic[3];
  PetscScalar ***g_vec;


  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
  assert(x[axis]==0 && n[axis]==N[axis]);

  DAVecGetArray(da, g, &g_vec);

  ic[0] = axis;
  ic[1] = (axis+1)%3;
  ic[2] = (axis+2)%3;

  i[ic[1]]=N[ic[1]]/2;
  i[ic[2]]=N[ic[2]]/2;
  for(i[ic[0]]=N[ic[0]]/2; i[ic[0]]<x[ic[0]]+n[ic[0]]; i[ic[0]]++)
    {
      fprintf(stderr,"%e\n",g_vec[i[2]][i[1]][i[0]]);
    }
  DAVecRestoreArray(da, g, &g_vec);
  fprintf(stderr,"end\n");
}

void Compute_g(BGY3dFourierData BDD, Vec g, Vec g0, Vec dg)
{
  DA da;
  int x[3], n[3], i[3];
  PetscScalar ***g_vec, ***dg_vec;


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

/*   VecView(g0,PETSC_VIEWER_STDERR_WORLD); */
/*   VecView(g,PETSC_VIEWER_STDERR_WORLD); */
}

void Compute_dg_kirk(BGY3dFourierData BDD, Vec g, Vec dg)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  FFT_DATA *(fg2_fft[3]), *g_fft, *dg_fft;
  real fac, k_fac;

  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];
  g_fft = BDD->g_fft;
  dg_fft = BDD->gfg2_fft;


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  VecCopy(g, dg);
  ShiftVec(da, dg, BDD->v[0], N);



  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], dg, BDD->f[dim]);
      ComputeFFTfromVec(da, BDD->fft_plan, BDD->v[dim], fg2_fft[dim]);
    }

  /* fft(g-1) */
  VecCopy(g, dg);
  VecShift(dg, -1.0);
  ComputeFFTfromVec(da, BDD->fft_plan, dg, g_fft);


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
	      dg_fft[index].re = 0;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k_fac = 1./(SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));


	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac *
		(fg2_fft[dim][index].re * g_fft[index].im
		 +fg2_fft[dim][index].im * g_fft[index].re) ;
	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac *
		(-fg2_fft[dim][index].re * g_fft[index].re
		 +fg2_fft[dim][index].im * g_fft[index].im);

	    }
	  index++;
	}

/*   for(i[0]=0;i[0]<n[0]*n[1]*n[2]; i[0]++) */
/*     fprintf(stderr,"%e\n", SQR(g_fft[i[0]].re)+SQR(g_fft[i[0]].im)); */
/*   exit(1); */

  ComputeVecfromFFT(da, BDD->fft_plan, dg, dg_fft);

  fac = -(PD->interval[1]-PD->interval[0])/(2.*M_PI);
  VecScale(dg, PD->h[0]*PD->h[1]*PD->h[2]*
	   PD->rho*PD->beta*fac
	   /PD->N[0]/PD->N[1]/PD->N[2]);



/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

void Compute_dg(BGY3dFourierData BDD, Vec g, Vec dg)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  FFT_DATA *(fg2_fft[3]), *g_fft, *dg_fft;
  real fac, k_fac;

  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];
  g_fft = BDD->g_fft;
  dg_fft = BDD->gfg2_fft;


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* fft(g-1) */
  VecCopy(g, dg);
  VecShift(dg, -1.0);
  ComputeFFTfromVec(da, BDD->fft_plan, dg, g_fft);


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
	      dg_fft[index].re = 0;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k_fac = 1./(SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));


	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac *
		(fg2_fft[dim][index].re * g_fft[index].im
		 +fg2_fft[dim][index].im * g_fft[index].re) ;
	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac *
		(-fg2_fft[dim][index].re * g_fft[index].re
		 +fg2_fft[dim][index].im * g_fft[index].im);

	    }
	  index++;
	}

/*   for(i[0]=0;i[0]<n[0]*n[1]*n[2]; i[0]++) */
/*     fprintf(stderr,"%e\n", SQR(g_fft[i[0]].re)+SQR(g_fft[i[0]].im)); */
/*   exit(1); */

  ComputeVecfromFFT(da, BDD->fft_plan, dg, dg_fft);

  fac = -(PD->interval[1]-PD->interval[0])/(2.*M_PI);
  VecScale(dg, PD->h[0]*PD->h[1]*PD->h[2]*
	   PD->rho*PD->beta*fac
	   /PD->N[0]/PD->N[1]/PD->N[2]);


/*   PetscScalar ***dg_vec; */
/*   DAVecGetArray(da, dg, &dg_vec); */
/*   PetscPrintf(PETSC_COMM_WORLD,"%e\n", dg_vec[63][63][63]); */
/*   DAVecRestoreArray(da, dg, &dg_vec); */


/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}





/* solve with product ansatz g=g0*dg */
Vec BGY3dDiv_solve_Fourier(ProblemData *PD, Vec g_ini, int vdim)
{
  BGY3dFourierData BDD;
  Vec g0, dg, dg_new, f, g;
  real a=0.9;
  int max_iter=25, iter, np;
  PetscScalar dg_norm, norm_tol=1.0e-6, norm;
  PetscTruth kflg;
  PetscViewer viewer;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dDiv equation with Fourier ansatz...\n");

  kflg = bgy3d_getopt_test ("-kirkwood");
  if(kflg)
    {
      PetscPrintf(PETSC_COMM_WORLD,"(Kirkwood)\n");
      MPI_Comm_size(PETSC_COMM_WORLD, &np);
      if(np>1)
	{
	  PetscPrintf(PETSC_COMM_WORLD,"Kirkwood only works with single process!\n");
	  return PETSC_NULL;
	}
      BDD = BGY3dFourierData_kirk_malloc(PD);
    }
  else
    {
      BDD = BGY3dFourierData_malloc(PD);
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


  VecDuplicate(BDD->g_ini, &g0);
  VecDuplicate(BDD->g_ini, &dg);
  VecDuplicate(BDD->g_ini, &dg_new);
  VecDuplicate(BDD->g_ini, &f);
  VecDuplicate(BDD->g_ini, &g);


  /* set initial guess*/
  VecSet(dg,0.0);
  VecSet(dg_new,0.0);
  VecCopy(BDD->g_ini, g0);


/*   VecView(v,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */


  for(iter=0; iter<max_iter; iter++)
    {


      /* g=g0+exp(-dg) */
      Compute_g(BDD, g, g0, dg);
/*       VecView(g,PETSC_VIEWER_STDERR_WORLD);   */
/*       VecView(g0,PETSC_VIEWER_STDERR_WORLD);   */
      /* f=integral(g) */
      if(kflg)
	Compute_dg_kirk(BDD, g, dg_new);
      else
	Compute_dg(BDD, g, dg_new);




      /* test for convergence of fixpoint iteration */
      VecWAXPY(f, -1.0, dg_new, dg);
      VecNorm(f, NORM_INFINITY, &dg_norm);
      VecNorm(f, NORM_2, &norm);



      PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norm: %e  %e\n", iter+1, dg_norm,
		  norm/PD->N[0]/PD->N[1]/PD->N[2]);

      if(dg_norm<norm_tol)
	{
	  VecAXPBY(dg, a, (1-a), dg_new);
	  break;
	}
      else
	/* simple mixing */
	VecAXPBY(dg, a, (1-a), dg_new);


    }

  /* g=g0*exp(-dg) */
  //Compute_g(BDD, g, g0, dg);
  VecCopy(dg, g);


/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
  if( kflg)
    {
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2_new.bin",
			    FILE_MODE_WRITE,&viewer);
      VecCopy(g,dg);
      ShiftVec(BDD->da, dg, BDD->v[0], PD->N);
      VecView(dg,viewer);
      PetscViewerDestroy(viewer);
      PetscPrintf(PETSC_COMM_WORLD,"g2 vector written to file \"g2_new.bin\".\n");
    }

  VecDestroy(g0);
  VecDestroy(dg);
  VecDestroy(dg_new);
  VecDestroy(f);

  // ExtractAxis(BDD, g, 0);


  BGY3dFourierData_free(BDD);


  return g;
}
