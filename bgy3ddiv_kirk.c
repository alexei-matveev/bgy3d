/*==========================================================*/
/*  $Id: bgy3ddiv.c,v 1.9 2006-08-21 16:22:15 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"

BGY3dDivData BGY3dDivData_kirk_malloc(ProblemData *PD, PetscTruth flg)
{
  BGY3dDivData BDD;
  DA da;
  real interval[2], h[3], N[3], L, r[3], r_s, **x_M, beta, h_g2;
  int i[3], x[3], n[3], dim, N_M, k, ic1, ic2, ic3, N_g2, index;
  PetscScalar ***ddU_vec, ***(f_vec[3]), ***boundary_vec, ***gini_vec, ***(v2_vec[3]), ***gSA_vec;
  PetscScalar *g2_vec;
  Vec g2;
  PetscViewer pview;
  int bufsize;

  BDD = (BGY3dDivData) malloc(sizeof(*BDD));


  BDD->LJ_params = (void* ) malloc(sizeof(real)*2);
  ((real*)(BDD->LJ_params))[0] = 1.0;   /* espilon */
  ((real*)(BDD->LJ_params))[1] = 1.0;   /* sigma   */

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
  PetscViewerBinaryOpen(PETSC_COMM_SELF,"g2file.bin" , PETSC_FILE_RDONLY,
			&pview);
  VecLoad( pview, VECSEQ, &g2);
  VecGetSize(g2,&N_g2);
  h_g2 = L/N_g2;



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
  VecGetArray(g2, &g2_vec);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */


	  FOR_DIM
	    {
	      r[dim] = i[dim]*h[dim];
	      if( i[dim]>=N[dim]/2)
		r[dim] -= L;
	    }

	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  FOR_DIM
	    f_vec[dim][i[2]][i[1]][i[0]] +=
	    Lennard_Jones_grad( r_s, r[dim], BDD->LJ_params);



	  ddU_vec[i[2]][i[1]][i[0]] +=
	    Lennard_Jones_ddU( r_s, r[0], BDD->LJ_params)+
	    Lennard_Jones_ddU( r_s, r[1], BDD->LJ_params)+
	    Lennard_Jones_ddU( r_s, r[2], BDD->LJ_params);

	  gini_vec[i[2]][i[1]][i[0]] *=
	    exp(-beta* Lennard_Jones( r_s, BDD->LJ_params));


	  index =(int) floor(r_s/h_g2);
	  if(index >= N_g2-1)
	    gSA_vec[i[2]][i[1]][i[0]] *= 1.0;
	  else
	    gSA_vec[i[2]][i[1]][i[0]] *=
	      g2_vec[index]+ (g2_vec[index+1]-g2_vec[index])/h_g2*
	      (r_s-index*h_g2);





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




/*   VecView(BDD->v2[0],PETSC_VIEWER_STDERR_WORLD);    */
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

void ComputeIntegralPart_kirk(BGY3dDivData BDD, Vec g, Vec f)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], dim, k, max_k;
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

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BDD->v[dim], g, BDD->f[dim]);
      ComputeFFTfromVec(da, BDD->fft_plan, BDD->v[dim], fg2_fft[dim], x, n, 0);
    }

  /* fft(g) */
  ComputeFFTfromVec(da, BDD->fft_plan, g, g_fft, x, n, 0);

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
      ComputeVecfromFFT(da, BDD->fft_plan, BDD->i[dim], gfg2_fft,
			x, n, 0.0);

      VecScale(BDD->i[dim], PD->h[0]*PD->h[1]*PD->h[2]*
	       PD->rho*PD->beta/PD->g_xm
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
/*   exit(1);   */

}

