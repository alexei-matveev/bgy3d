/*==========================================================*/
/*  $Id: bgy3ddiv.c,v 1.28 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fft.h"
#include "bgy3ddiv.h"


// real Lennard_Jones_ddU(real r, real xr, void *LJ_params)
real Lennard_Jones_ddU(real r, real xr, real epsilon, real sigma)
{
  // real epsilon, sigma, sr6, r2, sr, re;
  real sr6, r2, sr, re;

  r += SHIFT;

  if(r==0)
    return +CUTOFF;

  // epsilon = ((double*)LJ_params)[0];
  // sigma   = ((double*)LJ_params)[1];

  r2 = 1./r/r;
  sr = sigma/r;
  sr6 = SQR(sr)*SQR(sr)*SQR(sr);

  re = 24.*epsilon*sr6*r2*(28.*sr6-8.)*xr*xr*r2
    -  24.*epsilon*sr6*r2*(2.*sr6-1.);

  if(fabs(re)>CUTOFF)
    return +CUTOFF;
  else
    return re;
}


BGY3dDivData BGY3dDivData_malloc(ProblemData *PD, PetscTruth flg)
{
  BGY3dDivData BDD;
  DA da;
  real interval[2], h[3], N[3], L, r[3], r_s, **x_M, beta, h_g2;
  int i[3], x[3], n[3], N_M, k, ic1, ic2, ic3, N_g2, index;
  PetscScalar ***ddU_vec, ***(f_vec[3]), ***boundary_vec, ***gini_vec, ***(v2_vec[3]), ***gSA_vec;
  PetscScalar *g2_vec;
  Vec g2;
  PetscViewer pview;
  int bufsize;
  real epsilon, sigma;

  BDD = (BGY3dDivData) malloc(sizeof(*BDD));


  // BDD->LJ_params = (void* ) malloc(sizeof(real)*2);
  // ((real*)(BDD->LJ_params))[0] = 1.0;   /* espilon */
  // ((real*)(BDD->LJ_params))[1] = 1.0;   /* sigma   */
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
  DACreateGlobalVector(da, &(BDD->ddU));
  DACreateGlobalVector(da, &(BDD->boundary));
  DACreateGlobalVector(da, &(BDD->g_ini));
  DACreateGlobalVector(da, &(BDD->g_SA));
  FOR_DIM
    {
      VecDuplicate(BDD->ddU, &(BDD->f[dim]));
      VecDuplicate(BDD->ddU, &(BDD->v2[dim]));
      VecDuplicate(BDD->ddU, &(BDD->v[dim]));
      VecDuplicate(BDD->ddU, &(BDD->i[dim]));
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

  /* Load g2 from file, sequential only */
  /* g2 has to be on a grid [0,L] with L=interval[1]-interval[0] */
  PetscViewerBinaryOpen(PETSC_COMM_SELF,"g2file.bin" , FILE_MODE_READ,
			&pview);
  VecLoad( pview, VECSEQ, &g2);
  PetscViewerDestroy(pview);
  VecGetSize(g2,&N_g2);
  h_g2 = L/N_g2;
  VecGetArray(g2, &g2_vec);
  //g2_vec[0] = g2_vec[1]+g2_vec[1]-g2_vec[2];



  VecSet(BDD->ddU,0.0);
  FOR_DIM
    VecSet(BDD->f[dim],0.0);
  VecSet(BDD->boundary, 0.0);

  VecSet(BDD->g_ini, 1.0);
  VecSet(BDD->g_SA, 1.0);

  DAVecGetArray(da, BDD->ddU, &ddU_vec);
  DAVecGetArray(da, BDD->boundary, &boundary_vec);
  DAVecGetArray(da, BDD->g_ini, &gini_vec);
  DAVecGetArray(da, BDD->g_SA, &gSA_vec);
  FOR_DIM
    {
      DAVecGetArray(da, BDD->f[dim], &(f_vec[dim]));
      DAVecGetArray(da, BDD->v2[dim], &(v2_vec[dim]));
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
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim]+interval[0]-x_M[k][dim];
		}

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	      FOR_DIM
		  f_vec[dim][i[2]][i[1]][i[0]] +=
		    // Lennard_Jones_grad( r_s, r[dim], BDD->LJ_params);
		    Lennard_Jones_grad( r_s, r[dim], epsilon, sigma);


	      /* falsch !!! */
	      ddU_vec[i[2]][i[1]][i[0]] +=
		// Lennard_Jones_ddU( r_s, r[0], BDD->LJ_params)+
		// Lennard_Jones_ddU( r_s, r[1], BDD->LJ_params)+
		// Lennard_Jones_ddU( r_s, r[2], BDD->LJ_params);
		Lennard_Jones_ddU( r_s, r[0], epsilon, sigma)+
		Lennard_Jones_ddU( r_s, r[1], epsilon, sigma)+
		Lennard_Jones_ddU( r_s, r[2], epsilon, sigma);
	      gini_vec[i[2]][i[1]][i[0]] *=
		exp(-beta* Lennard_Jones( r_s, epsilon, sigma));


	      index =(int) floor(r_s/h_g2);
	      if(index >= N_g2-1)
		gSA_vec[i[2]][i[1]][i[0]] *= 1.0;
	      else
		gSA_vec[i[2]][i[1]][i[0]] *=
		  g2_vec[index]+ (g2_vec[index+1]-g2_vec[index])/h_g2*
		  (r_s-index*h_g2);

	    }




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
  DAVecRestoreArray(da, BDD->ddU, &ddU_vec);
  DAVecRestoreArray(da, BDD->boundary, &boundary_vec);
  DAVecRestoreArray(da, BDD->g_ini, &gini_vec);
  DAVecRestoreArray(da, BDD->g_SA, &gSA_vec);
  FOR_DIM
    {
      DAVecRestoreArray(da, BDD->f[dim], &(f_vec[dim]));
      DAVecRestoreArray(da, BDD->v2[dim], &(v2_vec[dim]));
    }




/*   VecView(BDD->g_SA,PETSC_VIEWER_STDERR_WORLD);     */
/*   exit(1);   */

  /* Create Matrix with appropriate non-zero structure */
  if(flg)
    DAGetMatrix( da, MATSEQAIJ, &(BDD->M));
  else
    DAGetMatrix( da, MATMPIAIJ, &(BDD->M));
  AssembleMatrix(BDD, da, BDD->M);

  FOR_DIM
    {
      DAGetMatrix( da, MATMPIAIJ, &(BDD->FD[dim]));
      AssembleFDMatrix(BDD, da, BDD->FD[dim], dim);
    }


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
					    BDD->v2[dim],
					    BDD->fg2_fft[dim]);
    }
  BDD->g_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));
  BDD->gfg2_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));


  /* set boundary vector for g0 in v2:*/
  FOR_DIM
    {
      VecSet(BDD->v2[dim],0.0);
      DAVecGetArray(da, BDD->v2[dim], &(v2_vec[dim]));
    }
  FOR_DIM
    {
      ic1=(dim)%3;
      ic2=(dim+1)%3;
      ic3=(dim+2)%3;
      if( x[ic1] == 0 )
	{
	  i[ic1] = 0;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      v2_vec[dim][i[2]][i[1]][i[0]] +=
		-0.5/h[dim];
	}
      if( x[ic1]+n[ic1] == PD->N[ic1])
	{
	  i[ic1] = PD->N[ic1]-1;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      v2_vec[dim][i[2]][i[1]][i[0]] +=
		0.5/h[dim];
	}
    }
  FOR_DIM
    DAVecRestoreArray(da, BDD->v2[dim], &(v2_vec[dim]));


  Molecule_free(x_M, N_M);
  VecDestroy(g2);

  return BDD;
}


BGY3dDivData BGY3dDivData_kirk_malloc(ProblemData *PD, PetscTruth flg)
{
  BGY3dDivData BDD;
  DA da;
  real interval[2], h[3], N[3], L, r[3], r_s, **x_M, beta, h_g2;
  int i[3], x[3], n[3], N_M, ic1, ic2, ic3, N_g2, index;
  PetscScalar ***ddU_vec, ***(f_vec[3]), ***boundary_vec, ***gini_vec, ***(v2_vec[3]), ***gSA_vec;
  PetscScalar *g2_vec;
  Vec g2;
  PetscViewer pview;
  int bufsize;
  real epsilon, sigma;

  BDD = (BGY3dDivData) malloc(sizeof(*BDD));


  // BDD->LJ_params = (void* ) malloc(sizeof(real)*2);
  // ((real*)(BDD->LJ_params))[0] = 1.0;   /* espilon */
  // ((real*)(BDD->LJ_params))[1] = 1.0;   /* sigma   */
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
  DACreateGlobalVector(da, &(BDD->ddU));
  DACreateGlobalVector(da, &(BDD->boundary));
  DACreateGlobalVector(da, &(BDD->g_ini));
  DACreateGlobalVector(da, &(BDD->g_SA));
  FOR_DIM
    {
      VecDuplicate(BDD->ddU, &(BDD->f[dim]));
      VecDuplicate(BDD->ddU, &(BDD->v2[dim]));
      VecDuplicate(BDD->ddU, &(BDD->v[dim]));
      VecDuplicate(BDD->ddU, &(BDD->i[dim]));
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

  /* Load g2 from file, sequential only */
  /* g2 has to be on a grid [0,L] with L=interval[1]-interval[0] */
  PetscViewerBinaryOpen(PETSC_COMM_SELF,"g2file.bin" , FILE_MODE_READ,
			&pview);
  VecLoad( pview, VECSEQ, &g2);
  PetscViewerDestroy(pview);
  VecGetSize(g2,&N_g2);
  h_g2 = L/N_g2;
  VecGetArray(g2, &g2_vec);
  //g2_vec[0] = g2_vec[1]+g2_vec[1]-g2_vec[2];


  VecSet(BDD->ddU,0.0);
  FOR_DIM
    VecSet(BDD->f[dim],0.0);
  VecSet(BDD->boundary, 0.0);

  VecSet(BDD->g_ini, 1.0);
  VecSet(BDD->g_SA, 1.0);

  DAVecGetArray(da, BDD->ddU, &ddU_vec);
  DAVecGetArray(da, BDD->boundary, &boundary_vec);
  DAVecGetArray(da, BDD->g_ini, &gini_vec);
  DAVecGetArray(da, BDD->g_SA, &gSA_vec);
  FOR_DIM
    {
      DAVecGetArray(da, BDD->f[dim], &(f_vec[dim]));
      DAVecGetArray(da, BDD->v2[dim], &(v2_vec[dim]));
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


	  index =(int) floor(r_s/h_g2);
	  if(index >= N_g2-1)
	    gSA_vec[i[2]][i[1]][i[0]] *= 1.0;
	  else
	    gSA_vec[i[2]][i[1]][i[0]] *=
	      g2_vec[index]+ (g2_vec[index+1]-g2_vec[index])/h_g2*
	      (r_s-index*h_g2);
	  ddU_vec[i[2]][i[1]][i[0]] =
	    exp(-beta*4./pow(r_s,6));

	  FOR_DIM
	    {
	      r[dim] = i[dim]*h[dim];
	      if( i[dim]>=N[dim]/2)
		r[dim] -= L;
	    }

	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  FOR_DIM
	    f_vec[dim][i[2]][i[1]][i[0]] +=
	    //Lennard_Jones_grad( r_s, r_s, BDD->LJ_params);
	    // Lennard_Jones_grad( r_s, r[dim], BDD->LJ_params);
	    Lennard_Jones_grad( r_s, r[dim], epsilon, sigma);



/* 	  ddU_vec[i[2]][i[1]][i[0]] +=  */
/* 	    Lennard_Jones_ddU( r_s, r[0], BDD->LJ_params)+ */
/* 	    Lennard_Jones_ddU( r_s, r[1], BDD->LJ_params)+ */
/* 	    Lennard_Jones_ddU( r_s, r[2], BDD->LJ_params); */
	}


  VecRestoreArray(g2, &g2_vec);
  DAVecRestoreArray(da, BDD->ddU, &ddU_vec);
  DAVecRestoreArray(da, BDD->boundary, &boundary_vec);
  DAVecRestoreArray(da, BDD->g_ini, &gini_vec);
  DAVecRestoreArray(da, BDD->g_SA, &gSA_vec);
  FOR_DIM
    {
      DAVecRestoreArray(da, BDD->f[dim], &(f_vec[dim]));
      DAVecRestoreArray(da, BDD->v2[dim], &(v2_vec[dim]));
    }



/*   VecView(BDD->f[0],PETSC_VIEWER_STDERR_WORLD);  */
/*   VecView(BDD->g_ini,PETSC_VIEWER_STDERR_WORLD);    */
/*   exit(1);  */

  /* Create Matrix with appropriate non-zero structure */
  if(flg)
    DAGetMatrix( da, MATSEQAIJ, &(BDD->M));
  else
    DAGetMatrix( da, MATMPIAIJ, &(BDD->M));
  AssembleMatrix(BDD, da, BDD->M);

  FOR_DIM
    {
      DAGetMatrix( da, MATMPIAIJ, &(BDD->FD[dim]));
      AssembleFDMatrix(BDD, da, BDD->FD[dim], dim);
    }


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
      BDD->fg2_fft[dim] = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));
    }

  BDD->g_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));
  BDD->gfg2_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));


  /* set boundary vector for g0 in v2:*/
  FOR_DIM
    {
      VecSet(BDD->v2[dim],0.0);
      DAVecGetArray(da, BDD->v2[dim], &(v2_vec[dim]));
    }
  FOR_DIM
    {
      ic1=(dim)%3;
      ic2=(dim+1)%3;
      ic3=(dim+2)%3;
      if( x[ic1] == 0 )
	{
	  i[ic1] = 0;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      v2_vec[dim][i[2]][i[1]][i[0]] +=
		-0.5/h[dim];
	}
      if( x[ic1]+n[ic1] == PD->N[ic1])
	{
	  i[ic1] = PD->N[ic1]-1;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      v2_vec[dim][i[2]][i[1]][i[0]] +=
		0.5/h[dim];
	}
    }
  FOR_DIM
    DAVecRestoreArray(da, BDD->v2[dim], &(v2_vec[dim]));


  Molecule_free(x_M, N_M);
  VecDestroy(g2);

  return BDD;
}



void BGY3dDivData_free(BGY3dDivData BDD)
{
  VecDestroy(BDD->ddU);
  FOR_DIM
    {
      VecDestroy(BDD->f[dim]);
      VecDestroy(BDD->v2[dim]);
      VecDestroy(BDD->v[dim]);
      VecDestroy(BDD->i[dim]);
      free(BDD->fg2_fft[dim]);
      MatDestroy(BDD->FD[dim]);
    }
  free(BDD->g_fft);
  free(BDD->gfg2_fft);
  VecDestroy(BDD->boundary);
  VecDestroy(BDD->g_ini);
  VecDestroy(BDD->g_SA);
  MatDestroy(BDD->M);
  DADestroy(BDD->da);
  // free(BDD->LJ_params);

  fft_3d_destroy_plan(BDD->fft_plan);

  free(BDD);
}



void AssembleMatrix(BGY3dDivData BDD, DA da, Mat M)
{
  ProblemData *PD;
  int x[3], n[3], i[3], N[3];
  MatStencil col[7],row;
  PetscScalar v[3], ***ddU_vec, ***(f_vec[3]);
  real h[3], beta;

  PD = BDD->PD;

  beta = BDD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  DAVecGetArray(da, BDD->ddU, &ddU_vec);
  FOR_DIM
    DAVecGetArray(da, BDD->f[dim], &(f_vec[dim]));

  MatZeroEntries(M);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  FOR_DIM
	    {
	      col[dim].i=i[0];
	      col[dim].j=i[1];
	      col[dim].k=i[2];
	    }
	  FOR_DIM
	    {
	      switch(dim)
		{
		case 0:
		  col[0].i -= 1;
		  col[2].i += 1;
		  break;
		case 1:
		  col[0].j -= 1;
		  col[2].j += 1;
		  break;
		case 2:
		  col[0].k -= 1;
		  col[2].k += 1;
		  break;
		}

	      /* values to enter */
	      v[0]= 1.0/SQR(h[dim]);
	      v[1]= -2.0/SQR(h[dim]);
	      v[2]= 1.0/SQR(h[dim]);
/* 	      v[0]= 1.0 ; */
/* 	      v[1]= -2.0 - SQR(h[dim])*SQR(beta)*SQR(f_vec[dim][i[2]][i[1]][i[0]]);  */
/* 	      v[2]= 1.0 ; */
	      /* check for boundary */
	      if( i[dim]== 0 )
		MatSetValuesStencil(M,1,&row,2,col+1,v+1,ADD_VALUES);
	      else if( i[dim]==N[dim]-1)
		MatSetValuesStencil(M,1,&row,2,col,v,ADD_VALUES);
	      else
		MatSetValuesStencil(M,1,&row,3,col,v,ADD_VALUES);

	      switch(dim)
		{
		case 0:
		  col[0].i = i[0];
		  col[2].i = i[0];
		  break;
		case 1:
		  col[0].j = i[1];
		  col[2].j = i[1];
		  break;
		case 2:
		  col[0].k = i[2];
		  col[2].k = i[2];
		  break;
		}
	    }
/* 	  v[0] = beta*ddU_vec[i[2]][i[1]][i[0]]; */
/* 	  MatSetValuesStencil(M,1,&row,1,col,v,ADD_VALUES); */


	}
  DAVecRestoreArray(da, BDD->ddU, &ddU_vec);
  FOR_DIM
    DAVecRestoreArray(da, BDD->f[dim], &(f_vec[dim]));

  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);

/*   MatView(M,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */

}




void AssembleFDMatrix(BGY3dDivData BDD, DA da, Mat M, int vdim)
{
  ProblemData *PD;
  int x[3], n[3], i[3], N[3];
  real h[3];
  MatStencil col[3],row;
  PetscScalar v[3];


  PD = BDD->PD;

  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  MatZeroEntries(M);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  FOR_DIM
	    {
	      col[dim].i=i[0];
	      col[dim].j=i[1];
	      col[dim].k=i[2];
	    }
	  switch(vdim)
	    {
	    case 0:
	      col[0].i -= 1;
	      col[2].i += 1;
	      break;
	    case 1:
	      col[0].j -= 1;
	      col[2].j += 1;
	      break;
	    case 2:
	      col[0].k -= 1;
	      col[2].k += 1;
	      break;
	    }

	  /* values to enter */
	  v[0]= -0.5/h[vdim] ;
	  v[1]= 0.0;
	  v[2]= 0.5/h[vdim] ;
	  /* check for boundary */
	  if( i[vdim]== 0 )
	    MatSetValuesStencil(M,1,&row,2,col+1,v+1,ADD_VALUES);
	  else if( i[vdim]==N[vdim]-1)
	    MatSetValuesStencil(M,1,&row,2,col,v,ADD_VALUES);
	  else
	    MatSetValuesStencil(M,1,&row,3,col,v,ADD_VALUES);


	}

  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);

/*   MatView(M,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */

}




void ComputeIntegralPart(BGY3dDivData BDD, Vec g, Vec f)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], k, max_k;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;


  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];
  g_fft = BDD->g_fft;
  gfg2_fft = BDD->gfg2_fft;


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
  max_k=n[0]*n[1]*n[2];

  /* fft(g) */
  ComputeFFTfromVec(da, BDD->fft_plan, g, BDD->g_fft);

  /* fft(g)*fft(fg2[dim]) */
  FOR_DIM
    {
      for(k=0; k<max_k; k++)
	{
	  gfg2_fft[k].re = fg2_fft[dim][k].re*g_fft[k].re
	    - fg2_fft[dim][k].im*g_fft[k].im;
	  gfg2_fft[k].im = fg2_fft[dim][k].re*g_fft[k].im
	    + fg2_fft[dim][k].im*g_fft[k].re;
	}

      /* i[dim]=fft^-1(fft(g)*fft(fg2)) */
      ComputeVecfromFFT(da, BDD->fft_plan, BDD->i[dim], gfg2_fft);

      VecScale(BDD->i[dim], PD->h[0]*PD->h[1]*PD->h[2]*
	       PD->rho*PD->beta
	       /PD->N[0]/PD->N[1]/PD->N[2]);
    }

  FOR_DIM
    {
      /* v= grad(Integral) */
      MatMult(BDD->FD[dim], BDD->i[dim], BDD->v[dim]);
    }

  /* f= v[0]+v[1]+v[2] */
  VecWAXPY(f, 1.0, BDD->v[0], BDD->v[1]);
  VecAXPY(f, 1.0, BDD->v[2]);

/*   VecView(g,PETSC_VIEWER_STDERR_WORLD);  */
/*   VecView(BDD->i[0],PETSC_VIEWER_STDERR_WORLD); */
/*   VecView(BDD->v[0],PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);   */

}

void ComputeIntegralPart_kirk(BGY3dDivData BDD, Vec g, Vec f)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], k, max_k;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;


  PD=BDD->PD;

  da = BDD->da;
  FOR_DIM
    fg2_fft[dim] = BDD->fg2_fft[dim];
  g_fft = BDD->g_fft;
  gfg2_fft = BDD->gfg2_fft;


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
  max_k=n[0]*n[1]*n[2];

  VecCopy(g, f);
  ShiftVec(da, f, BDD->v[0], PD->N);
  // ShiftVec(BDD, f);



  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], f, BDD->f[dim]);
      ComputeFFTfromVec(da, BDD->fft_plan, BDD->v[dim], fg2_fft[dim]);
    }

  /* fft(g) */
  VecCopy(g, f);
  VecShift(f, -1.0);
  ComputeFFTfromVec(da, BDD->fft_plan, f, g_fft);
/*   for(k=0; k<max_k; k++) */
/*     fprintf(stderr,"%e\n",g_fft[k].re); */
/*   ComputeVecfromFFT(da, BDD->fft_plan, f, g_fft,  */
/* 			x, n, 0.0); */
/*   VecScale(f,1./PD->N[0]/PD->N[1]/PD->N[2]); */
/*   VecCopy(BDD->v[0], f); */
/*   ShiftVec(BDD, f); */
/*   VecView(g,PETSC_VIEWER_STDERR_WORLD);   */
/*   VecView(BDD->v[0],PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */


  /* fft(g)*fft(fg2[dim]) */
  FOR_DIM
    {
      for(k=0; k<max_k; k++)
	{
	  gfg2_fft[k].re = fg2_fft[dim][k].re*g_fft[k].re
	    - fg2_fft[dim][k].im*g_fft[k].im;

	  gfg2_fft[k].im = fg2_fft[dim][k].re*g_fft[k].im
	    + fg2_fft[dim][k].im*g_fft[k].re;
	}

      /* i[dim]=fft^-1(fft(g)*fft(fg2)) */
      ComputeVecfromFFT(da, BDD->fft_plan, BDD->i[dim], gfg2_fft);

      VecScale(BDD->i[dim], PD->h[0]*PD->h[1]*PD->h[2]*
	       PD->rho*PD->beta/PD->g_xm
	       /PD->N[0]/PD->N[1]/PD->N[2]);
    }

  FOR_DIM
    {
      //VecPointwiseMult(BDD->i[dim], BDD->i[dim], BDD->ddU);
      /* v= grad(Integral) */
      MatMult(BDD->FD[dim], BDD->i[dim], BDD->v[dim]);

    }

  /* f= v[0]+v[1]+v[2] */
  VecWAXPY(f, 1.0, BDD->v[0], BDD->v[1]);
  VecAXPY(f, 1.0, BDD->v[2]);
/*   VecWAXPY(f, 1.0, BDD->i[0], BDD->i[1]); */
/*   VecAXPY(f, 1.0, BDD->i[2]); */

/*   VecView(BDD->i[0],PETSC_VIEWER_STDERR_WORLD);  */
/*   VecView(f,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);    */

}



void AssembleSystemMatrix(BGY3dDivData BDD, Mat SM, Vec f)
{
  DA da;
  ProblemData *PD;
  int x[3], n[3], i[3], N[3];
  real h[3];
  MatStencil col,row;
  PetscScalar ***f_vec;


  PD = BDD->PD;
  da = BDD->da;

  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  MatZeroEntries(SM);
  MatCopy(BDD->M, SM, SAME_NONZERO_PATTERN);

  DAVecGetArray(da, f, &f_vec);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  col.i=i[0];
	  col.j=i[1];
	  col.k=i[2];

	  MatSetValuesStencil(SM,1,&row,1,&col,&(f_vec[i[2]][i[1]][i[0]]),
			      ADD_VALUES);


	}
  DAVecRestoreArray(da, f, &f_vec);
  MatAssemblyBegin(SM,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(SM,MAT_FINAL_ASSEMBLY);

/*   MatView(M,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */

}

void AssembleSystemMatrix2(BGY3dDivData BDD, Mat SM)
{
  DA da;
  ProblemData *PD;
  int x[3], n[3], i[3], N[3];
  real h[3];
  MatStencil col[2],row;
  PetscScalar ***(i_vec[3]), val[2];


  PD = BDD->PD;
  da = BDD->da;

  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  MatZeroEntries(SM);
  MatCopy(BDD->M, SM, SAME_NONZERO_PATTERN);

  FOR_DIM
    DAVecGetArray(da, BDD->i[dim], &(i_vec[dim]));
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  FOR_DIM
	    {
	      col[0].i=i[0];
	      col[0].j=i[1];
	      col[0].k=i[2];
	      col[1].i=i[0];
	      col[1].j=i[1];
	      col[1].k=i[2];
	      switch(dim)
		{
		case 0:
		  col[0].i -= 1; col[1].i += 1;
		  if(i[0]==0)
		    val[1] =  0.5*i_vec[0][i[2]][i[1]][i[0]+1]/h[0];
		  else if(i[0] == N[0]-1)
		    val[0] = -0.5*i_vec[0][i[2]][i[1]][i[0]-1]/h[0];
		  else
		    {
		      val[0] = -0.5*i_vec[0][i[2]][i[1]][i[0]-1]/h[0];
		      val[1] =  0.5*i_vec[0][i[2]][i[1]][i[0]+1]/h[0];
		    }
		  break;
		case 1:
		  col[0].j -= 1; col[1].j += 1;
		  if(i[1]==0)
		    val[1] =  0.5*i_vec[1][i[2]][i[1]+1][i[0]]/h[0];
		  else if(i[1]==N[1]-1)
		    val[0] = -0.5*i_vec[1][i[2]][i[1]-1][i[0]]/h[0];
		  else
		    {
		      val[0] = -0.5*i_vec[1][i[2]][i[1]-1][i[0]]/h[0];
		      val[1] =  0.5*i_vec[1][i[2]][i[1]+1][i[0]]/h[0];
		    }
		  break;
		case 2:
		  col[0].k -= 1; col[1].k += 1;
		  if(i[2]==0)
		    val[1] =  0.5*i_vec[2][i[2]+1][i[1]][i[0]]/h[0];
		  else if(i[2]==N[2]-1)
		    val[0] = -0.5*i_vec[2][i[2]-1][i[1]][i[0]]/h[0];
		  else
		    {
		      val[0] = -0.5*i_vec[2][i[2]-1][i[1]][i[0]]/h[0];
		      val[1] =  0.5*i_vec[2][i[2]+1][i[1]][i[0]]/h[0];
		    }
		  break;
		}
	      if(i[dim]==0)
		MatSetValuesStencil(SM,1,&row,1,col+1,val+1, ADD_VALUES);
	      else if(i[dim]==N[dim]-1)
		MatSetValuesStencil(SM,1,&row,1,col,val, ADD_VALUES);
	      else
		MatSetValuesStencil(SM,1,&row,2,col,val, ADD_VALUES);
	    }

	}
  FOR_DIM
    DAVecRestoreArray(da, BDD->i[dim], &(i_vec[dim]));
  MatAssemblyBegin(SM,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(SM,MAT_FINAL_ASSEMBLY);

/*   MatView(M,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */

}



void AssembleSystemMatrix_part2(BGY3dDivData BDD, Mat SM)
{
  DA da;
  ProblemData *PD;
  int x[3], n[3], i[3], N[3];
  real h[3];
  MatStencil col[3],row;
  PetscScalar ***(v_vec[3]), val[3], ***ddU_vec;
  real beta;

  PD = BDD->PD;
  da = BDD->da;

  beta = PD->beta;

  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  FOR_DIM
    /* v=beta*grad(U)+integral */
    VecWAXPY(BDD->v[dim], PD->beta, BDD->f[dim], BDD->i[dim]);

  FOR_DIM
    DAVecGetArray(da, BDD->v[dim], &(v_vec[dim]));
  DAVecGetArray(da, BDD->ddU, &(ddU_vec));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  FOR_DIM
	    {
	      col[0].i=i[0];
	      col[0].j=i[1];
	      col[0].k=i[2];
	      col[1].i=i[0];
	      col[1].j=i[1];
	      col[1].k=i[2];
	      col[2].i=i[0];
	      col[2].j=i[1];
	      col[2].k=i[2];
	      switch(dim)
		{
		case 0:
		  col[0].i -= 1;
		  col[2].i += 1;
		  break;
		case 1:
		  col[0].j -= 1;
		  col[2].j += 1;
		  break;
		case 2:
		  col[0].k -= 1;
		  col[2].k += 1;
		  break;
		}
	      val[0] = -0.5/h[dim]*v_vec[dim][i[2]][i[1]][i[0]];
	      val[1] = beta*ddU_vec[i[2]][i[1]][i[0]];
	      val[2] = +0.5/h[dim]*v_vec[dim][i[2]][i[1]][i[0]];
	      /* check for boundary */
	      if( i[dim]== 0 )
		MatSetValuesStencil(SM,1,&row,2,col+1,val+1,ADD_VALUES);
	      else if( i[dim]==N[dim]-1)
		MatSetValuesStencil(SM,1,&row,2,col,val,ADD_VALUES);
	      else
		MatSetValuesStencil(SM,1,&row,3,col,val,ADD_VALUES);
	    }
	}
  FOR_DIM
    DAVecRestoreArray(da, BDD->v[dim], &(v_vec[dim]));
  DAVecRestoreArray(da, BDD->ddU, &(ddU_vec));
  MatAssemblyBegin(SM,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(SM,MAT_FINAL_ASSEMBLY);

/*   MatView(M,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */

}

void AssembleSystemMatrix_part2b(BGY3dDivData BDD, Mat SM)
{
  DA da;
  ProblemData *PD;
  int x[3], n[3], i[3], N[3];
  real h[3];
  MatStencil col[2],row;
  PetscScalar ***(i_vec[3]), val[2];


  PD = BDD->PD;
  da = BDD->da;


  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));




  FOR_DIM
    DAVecGetArray(da, BDD->i[dim], &(i_vec[dim]));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  FOR_DIM
	    {
	      col[0].i=i[0];
	      col[0].j=i[1];
	      col[0].k=i[2];
	      col[1].i=i[0];
	      col[1].j=i[1];
	      col[1].k=i[2];
	      switch(dim)
		{
		case 0:
		  col[0].i -= 1;
		  col[1].i += 1;
		  break;
		case 1:
		  col[0].j -= 1;
		  col[1].j += 1;
		  break;
		case 2:
		  col[0].k -= 1;
		  col[1].k += 1;
		  break;
		}
	      val[0] = -0.5/h[dim]*i_vec[dim][i[2]][i[1]][i[0]];
	      val[1] = +0.5/h[dim]*i_vec[dim][i[2]][i[1]][i[0]];
	      /* check for boundary */
	      if( i[dim]== 0 )
		MatSetValuesStencil(SM,1,&row,1,col+1,val+1,ADD_VALUES);
	      else if( i[dim]==N[dim]-1)
		MatSetValuesStencil(SM,1,&row,1,col,val,ADD_VALUES);
	      else
		MatSetValuesStencil(SM,1,&row,2,col,val,ADD_VALUES);
	    }
	}
  FOR_DIM
    DAVecRestoreArray(da, BDD->i[dim], &(i_vec[dim]));
  MatAssemblyBegin(SM,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(SM,MAT_FINAL_ASSEMBLY);

/*   MatView(M,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */

}


void ComputeRHS(BGY3dDivData BDD, Vec b, Vec g0, Vec f)
{
  /* v= grad(g)*integral  */
  FOR_DIM
    {
      MatMult(BDD->FD[dim], g0, BDD->v[dim]);
      VecAXPY(BDD->v[dim], 1.0, BDD->v2[dim]); /* boundary */
      VecPointwiseMult(BDD->v[dim], BDD->v[dim], BDD->i[dim]);
    }

  /* b= grad(g)'*integral */
  VecWAXPY(b, 1.0, BDD->v[0], BDD->v[1]);
  VecAXPY(b, 1.0, BDD->v[2]);



  /* b += grad(integral)*g0 */
  VecPointwiseMult(BDD->v[0], g0, f);
  VecAXPY(b, 1.0, BDD->v[0]);

  VecScale(b,-1.0);

/*   VecView(b,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);   */

}

void ComputeRHS2(BGY3dDivData BDD, Vec b)
{
  DA da;
  ProblemData *PD;
  int x[3], n[3], i[3], N[3], ic1, ic2, ic3;
  real h[3];
  PetscScalar ***(i_vec[3]), ***b_vec;


  PD = BDD->PD;
  da = BDD->da;


  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  VecSet(b,0.0);

  FOR_DIM
    DAVecGetArray(da, BDD->i[dim], &(i_vec[dim]));
  DAVecGetArray(da, b, &(b_vec));

  /* set boundary vector */
  FOR_DIM
    {
      ic1=(dim)%3;
      ic2=(dim+1)%3;
      ic3=(dim+2)%3;
      if( x[ic1] == 0 )
	{
	  i[ic1] = 0;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      b_vec[i[2]][i[1]][i[0]] +=
		-(1.0/SQR(h[dim])-0.5/h[dim]*i_vec[dim][i[2]][i[1]][i[0]]);
	}
      if( x[ic1]+n[ic1] == PD->N[ic1])
	{
	  i[ic1] = PD->N[ic1]-1;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      b_vec[i[2]][i[1]][i[0]] +=
		-(1.0/SQR(h[dim])+0.5/h[dim]*i_vec[dim][i[2]][i[1]][i[0]]);
	}
    }

  FOR_DIM
    DAVecRestoreArray(da, BDD->i[dim], &(i_vec[dim]));
  DAVecRestoreArray(da, b, &(b_vec));


}

void ComputeRHS_laplace(BGY3dDivData BDD, Vec b)
{
  DA da;
  ProblemData *PD;
  int x[3], n[3], i[3], N[3], ic1, ic2, ic3;
  real h[3];
  PetscScalar ***(i_vec[3]), ***b_vec;


  PD = BDD->PD;
  da = BDD->da;


  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  VecSet(b,0.0);

  FOR_DIM
    DAVecGetArray(da, BDD->i[dim], &(i_vec[dim]));
  DAVecGetArray(da, b, &(b_vec));

  /* set boundary vector */
  FOR_DIM
    {
      ic1=(dim)%3;
      ic2=(dim+1)%3;
      ic3=(dim+2)%3;
      if( x[ic1] == 0 )
	{
	  i[ic1] = 0;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      b_vec[i[2]][i[1]][i[0]] +=
		-(1.0/SQR(h[dim]));
	}
      if( x[ic1]+n[ic1] == PD->N[ic1])
	{
	  i[ic1] = PD->N[ic1]-1;
	  for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	    for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	      b_vec[i[2]][i[1]][i[0]] +=
		-(1.0/SQR(h[dim]));
	}
    }

  FOR_DIM
    DAVecRestoreArray(da, BDD->i[dim], &(i_vec[dim]));
  DAVecRestoreArray(da, b, &(b_vec));


}




void ShiftVec(DA da, Vec g, Vec scratch, int N[3])
{
  int x[3], n[3], i[3], ic[3], add[3];
  PetscScalar  ***g_vec, ***temp_vec;


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  VecCopy(g, scratch);

  FOR_DIM
    if( N[dim]%2)
      add[dim] = N[dim]/2+1;
    else
      add[dim] = N[dim]/2;

  DAVecGetArray(da, g, &(g_vec));
  DAVecGetArray(da, scratch, &(temp_vec));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  FOR_DIM
	    {
	      if(i[dim]<N[dim]/2)
		ic[dim] = i[dim] + add[dim];
	      else
		ic[dim] = i[dim] - N[dim]/2;
	    }

	  g_vec[ic[2]][ic[1]][ic[0]] = temp_vec[i[2]][i[1]][i[0]];
	}

  DAVecRestoreArray(da, g, &(g_vec));
  DAVecRestoreArray(da, scratch, &(temp_vec));
}

void ComputeBGY3dDiv_F(BGY3dDivData BDD, Mat SM, Vec g0, Vec dg, Vec g,
		       Vec b, Vec f)
{

  VecWAXPY(g, 1.0, dg, g0);
  /* f=integral(g) */
  ComputeIntegralPart(BDD, g, f);

  /* Assemble matrix to solve */
  AssembleSystemMatrix(BDD, SM, f);
  AssembleSystemMatrix_part2(BDD, SM);

  /* set right hand site */
  ComputeRHS(BDD, b, g0, f);

  MatMult(SM, dg, f);

  VecAXPY(f,-1.0, b);

/*   VecView(f,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */

}

void ComputeBGY3dDiv_F2(BGY3dDivData BDD, Mat SM, Vec g0, Vec dg, Vec g,
			Vec b, Vec f, PetscTruth kflg)
{

  VecPointwiseMult(g, dg, g0);
  /* f=integral(g) */
  if(kflg)
    ComputeIntegralPart_kirk(BDD, g, f);
  else
    ComputeIntegralPart(BDD, g, f);

  /* Assemble matrix to solve */
  AssembleSystemMatrix(BDD, SM, f);
  AssembleSystemMatrix_part2b(BDD, SM);

  /* set right hand site */
  if(kflg)
    ComputeRHS2(BDD, b);
  else
    ComputeRHS2(BDD, b);

  MatMult(SM, dg, f);

  VecAXPY(f,-1.0, b);

/*   VecView(f,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */

}

void ComputeBGY3dDiv_F3(BGY3dDivData BDD, Mat SM, Vec g0, Vec dg, Vec g,
			Vec b, Vec f, PetscTruth kflg)
{
  Vec help[3], part[3]; // was declared as part[dim];


  FOR_DIM
    VecDuplicate(g, &(help[dim]));
  FOR_DIM
    VecDuplicate(g, &(part[dim]));

  VecPointwiseMult(g, dg, g0);
  /* f=integral(g) */
  if(kflg)
    ComputeIntegralPart_kirk(BDD, g, f);
  else
    ComputeIntegralPart(BDD, g, f);

  MatZeroEntries(SM);
  MatCopy(BDD->M, SM, SAME_NONZERO_PATTERN);

  FOR_DIM
    MatMult(BDD->FD[dim], BDD->i[dim], BDD->v[dim]);

  VecWAXPY(part[1], 1.0, BDD->v[0], BDD->v[1]);
  VecAXPY(part[1], 1.0, BDD->v[2]);
  VecPointwiseMult(part[1], part[1], dg);

/*   VecPointwiseMult(BDD->v[0], BDD->v[0], BDD->g_SA); */
/*   VecView(BDD->g_SA,PETSC_VIEWER_STDERR_WORLD); */
/*   VecView(BDD->v[0],PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

  FOR_DIM
    {
      MatMult(BDD->FD[dim], dg, BDD->v[dim]);
      VecAXPY(BDD->v[dim], 1.0, BDD->v2[dim]);
      VecPointwiseMult(help[dim], BDD->v[dim], BDD->i[dim]);
    }
  VecWAXPY(part[2], 1.0, help[0], help[1]);
  VecAXPY(part[2], 1.0, help[2]);

  MatMult(SM, dg, part[0]);
  ComputeRHS_laplace(BDD, b);
  VecAXPY(part[0],-1.0, b);

  VecWAXPY(f, 1.0, part[0], part[1]);
  VecAXPY(f, 1.0, part[2]);


/*   VecView(f,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */
  VecView(part[0],PETSC_VIEWER_STDERR_WORLD);
  VecView(part[1],PETSC_VIEWER_STDERR_WORLD);
  VecView(part[2],PETSC_VIEWER_STDERR_WORLD);

  FOR_DIM
    VecDestroy(help[dim]);
  FOR_DIM
    VecDestroy(part[dim]);

}


void ComputeBGY3d_F(BGY3dDivData BDD, Vec g0, Vec dg, Vec g,
		    Vec b, Vec f, PetscTruth kflg, int vdim)
{

  VecPointwiseMult(g, dg, g0);


  /* f=integral(g) */
  if(kflg)
    ComputeIntegralPart_kirk(BDD, g, f);
  else
    ComputeIntegralPart(BDD, g, f);


  /* set right hand site */
  VecCopy(BDD->v2[vdim], b);
  VecScale(b,-1.0);





/*   VecCopy(g,f); */
/*   ShiftVec(BDD, f); */
/*   VecPointwiseMult(BDD->v[2], f, BDD->f[vdim]); */
/*   VecScale(BDD->v[2], BDD->beta); */
/*   ShiftVec(BDD,BDD->v[2]); */

  VecPointwiseMult(BDD->v[0], dg, BDD->i[vdim]);
  MatMult(BDD->FD[vdim], dg, BDD->v[1]);


  VecWAXPY(f, 1.0, BDD->v[0], BDD->v[1]);
/*   VecAXPY(f, 1.0, BDD->v[2]);  */
  VecAXPY(f,-1.0, b);

/*   VecAXPY(BDD->v[1],-1.0, b); */
/*   VecView(BDD->v[0],PETSC_VIEWER_STDERR_WORLD);   */
/*   VecView(BDD->v[1],PETSC_VIEWER_STDERR_WORLD);   */
/*   VecView(BDD->i[0],PETSC_VIEWER_STDERR_WORLD);   */
  //exit(1);

}


Vec BGY3dDiv_solve(ProblemData *PD, Vec g_ini, int vdim)
{
  BGY3dDivData BDD;
  Vec g0, dg, dg_new, b, f, g;
  Mat SM;
  KSP ksp;
  PC pc;
  real a=0.5;
  int max_iter=25, iter, l_iter;
  PetscScalar dg_norm, f_norm, r_norm, norm_tol=1.0e-6;
  PetscTruth flg;


  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dDiv equation...\n");

  flg = bgy3d_getopt_test ("-seq");

  BDD = BGY3dDivData_malloc(PD, flg);

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
  if(flg)
    PetscPrintf(PETSC_COMM_WORLD,"Using sequential matrix format.\n");
  PetscPrintf(PETSC_COMM_WORLD,"lambda = %f\n",a);
  PetscPrintf(PETSC_COMM_WORLD,"tolerance = %e\n",norm_tol);
  PetscPrintf(PETSC_COMM_WORLD,"max_iter = %d\n",max_iter);


  VecDuplicate(BDD->ddU, &g0);
  VecDuplicate(BDD->ddU, &b);
  VecDuplicate(BDD->ddU, &dg);
  VecDuplicate(BDD->ddU, &dg_new);
  VecDuplicate(BDD->ddU, &f);
  VecDuplicate(BDD->ddU, &g);
  if(flg)
    DAGetMatrix( BDD->da, MATSEQAIJ, &SM);
  else
    DAGetMatrix( BDD->da, MATMPIAIJ, &SM);
  MatSetOption(SM, MAT_NO_NEW_NONZERO_LOCATIONS);

  /* set initial guess*/
  VecSet(dg,0.0);
  VecSet(dg_new,0.0);
  VecCopy(BDD->g_ini, g0);


/*   VecView(v,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */
  /**********************************/
  /* initialize ksp */
  /**********************************/
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPGetPC(ksp, &pc);
  PCSetType( pc, PCJACOBI);

  KSPSetType(ksp, KSPGMRES);

  /* set rtol, atol, dtol, maxits */
  // KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  KSPSetTolerances(ksp, 1.0e-6, 1.0e-50, 1.0e+5, 1000);

  /* runtime options will override options set above */
  KSPSetFromOptions(ksp);
  /**********************************/


  for(iter=0; iter<max_iter; iter++)
    {

      VecSet(g,0.0);
      VecSet(b,0.0);
      VecSet(f,0.0);
      VecSet(dg_new,0.0);

      /* g=g0+dg */
      VecWAXPY(g, 1.0, dg, g0);
      /* f=integral(g) */
      ComputeIntegralPart(BDD, g, f);

      /* Assemble matrix to solve */
      AssembleSystemMatrix(BDD, SM, f);
      AssembleSystemMatrix_part2(BDD, SM);

      /* set right hand site */
      ComputeRHS(BDD, b, g0, f);


      /* set matrix in ksp */
      KSPSetOperators(ksp, SM, SM, SAME_NONZERO_PATTERN);


      /* solve system */
      VecSet(dg_new, 0.0);
      KSPSolve(ksp, b, dg_new);

      KSPGetResidualNorm(ksp, &r_norm);
      KSPGetIterationNumber(ksp, &l_iter);

/*       MatView(SM,PETSC_VIEWER_STDERR_WORLD); */
/*       exit(1); */

      /* test for convergence of fixpoint iteration */
      VecWAXPY(b, -1.0, dg_new, dg);
      VecNorm(b, NORM_2, &dg_norm);

      ComputeBGY3dDiv_F(BDD, SM, g0, dg_new, g, b, f);
      VecNorm(f, NORM_2, &f_norm);


      PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg, f function norm: %e, %e\tRnorm= %e after %d iterations.\n",
		  iter+1, dg_norm, f_norm, r_norm, l_iter);

      if(dg_norm<norm_tol)
	{
	  VecAXPBY(dg, a, (1-a), dg_new);
	  break;
	}
      else
	/* simple mixing */
	VecAXPBY(dg, a, (1-a), dg_new);
    }

  /* g=g0+dg */
  VecWAXPY(g, 1.0, dg, g0);
  //VecCopy(dg, g);

  VecDestroy(g0);
  VecDestroy(dg);
  VecDestroy(dg_new);
  VecDestroy(b);
  VecDestroy(f);
  MatDestroy(SM);

  KSPDestroy(ksp);
  BGY3dDivData_free(BDD);

  return g;
}



/* solve with product ansatz g=g0*dg */
Vec BGY3dDiv_solve2(ProblemData *PD, Vec g_ini, int vdim)
{
  BGY3dDivData BDD;
  Vec g0, dg, dg_new, b, f, g;
  Mat SM;
  KSP ksp;
  PC pc;
  real a=0.9;
  int max_iter=25, iter, l_iter, np;
  PetscScalar dg_norm, f_norm, r_norm, norm_tol=1.0e-6; // f2_norm
  PetscTruth flg, kflg;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dDiv equation with product ansatz...");

  flg = bgy3d_getopt_test ("-seq");
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
      BDD = BGY3dDivData_kirk_malloc(PD,flg);
    }
  else
    {
      BDD = BGY3dDivData_malloc(PD,flg);
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
  if(flg)
    PetscPrintf(PETSC_COMM_WORLD,"Using sequential matrix format.\n");
  PetscPrintf(PETSC_COMM_WORLD,"lambda = %f\n",a);
  PetscPrintf(PETSC_COMM_WORLD,"tolerance = %e\n",norm_tol);
  PetscPrintf(PETSC_COMM_WORLD,"max_iter = %d\n",max_iter);


  VecDuplicate(BDD->ddU, &g0);
  VecDuplicate(BDD->ddU, &b);
  VecDuplicate(BDD->ddU, &dg);
  VecDuplicate(BDD->ddU, &dg_new);
  VecDuplicate(BDD->ddU, &f);
  VecDuplicate(BDD->ddU, &g);
  if(flg)
    DAGetMatrix( BDD->da, MATSEQAIJ, &SM);
  else
    DAGetMatrix( BDD->da, MATMPIAIJ, &SM);
  MatSetOption(SM, MAT_NO_NEW_NONZERO_LOCATIONS);

  /* set initial guess*/
  VecSet(dg,1.0);
  VecSet(dg_new,1.0);
  VecCopy(BDD->g_ini, g0);


/*   VecView(v,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */
  /**********************************/
  /* initialize ksp */
  /**********************************/
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPGetPC(ksp, &pc);
  PCSetType( pc, PCJACOBI);

  KSPSetType(ksp, KSPBCGS);

  /* set rtol, atol, dtol, maxits */
  // KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  KSPSetTolerances(ksp, 1.0e-10, 1.0e-50, 1.0e+5, 1000);

  /* runtime options will override options set above */
  KSPSetFromOptions(ksp);
  /**********************************/


  for(iter=0; iter<max_iter; iter++)
    {


      /* g=g0+dg */
      VecPointwiseMult(g, dg, g0);
      /* f=integral(g) */
      if(kflg)
	ComputeIntegralPart_kirk(BDD, g, f);
      else
	ComputeIntegralPart(BDD, g, f);

      /* Assemble matrix to solve */
      AssembleSystemMatrix(BDD, SM, f);
      AssembleSystemMatrix_part2b(BDD, SM);

      //AssembleSystemMatrix2(BDD, SM);

/*       MatView(SM,PETSC_VIEWER_STDERR_WORLD); */
/*       exit(1); */

      /* set right hand site */
      ComputeRHS2(BDD, b);
      //ComputeRHS_laplace(BDD, b);

      /* set matrix in ksp */
      KSPSetOperators(ksp, SM, SM, SAME_NONZERO_PATTERN);


      /* solve system */
      //VecSet(dg_new, 1.0);
      KSPSolve(ksp, b, dg_new);

      KSPGetResidualNorm(ksp, &r_norm);
      KSPGetIterationNumber(ksp, &l_iter);



      /* test for convergence of fixpoint iteration */
      VecWAXPY(b, -1.0, dg_new, dg);
      VecNorm(b, NORM_2, &dg_norm);

      ComputeBGY3dDiv_F2(BDD, SM, g0, dg_new, g, b, f, kflg);
      VecNorm(f, NORM_2, &f_norm);

      /**************************/

      //ComputeBGY3d_F(BDD, g0, dg, g, b, f, kflg, 0);
      //ComputeBGY3dDiv_F3(BDD, SM, g0, dg_new, g, b, f, kflg);
      //VecNorm(f, NORM_2, &f2_norm);


      PetscPrintf(PETSC_COMM_WORLD, "iter %d: dg, f function norm: %e, %e\tRnorm= %e after %d iterations.\n",
		  iter+1, dg_norm, f_norm, r_norm, l_iter);

      if(f_norm<norm_tol)
	{
	  VecAXPBY(dg, a, (1-a), dg_new);
	  break;
	}
      else
	/* simple mixing */
	VecAXPBY(dg, a, (1-a), dg_new);


    }

  /* g=g0+dg */
  VecPointwiseMult(g, dg, g0);
  //VecCopy(dg, g);
/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
/*   VecWAXPY(b, -1.0, g, BDD->g_SA); */
/*   VecNorm(b, NORM_2, &dg_norm); */
/*   PetscPrintf(PETSC_COMM_WORLD,"Difference to SA: f norm= %e\n",dg_norm);  */

/*   ComputeBGY3d_F(BDD, g0, dg, g, b, f, kflg, 0); */
  //ComputeBGY3dDiv_F3(BDD, SM, g0, BDD->g_SA, g, b, f, kflg);
/*   VecNorm(f, NORM_INFINITY, &f2_norm); */
/*   PetscPrintf(PETSC_COMM_WORLD,"norm= %e\n", f2_norm); */
  //VecPointwiseMult(g, BDD->g_SA, g0);
  //VecCopy(BDD->g_SA, g);
/*   if(kflg) */
/*     ShiftVec(BDD,g); */

  VecDestroy(g0);
  VecDestroy(dg);
  VecDestroy(dg_new);
  VecDestroy(b);
  VecDestroy(f);
  MatDestroy(SM);

  KSPDestroy(ksp);
  BGY3dDivData_free(BDD);



  return g;
}
