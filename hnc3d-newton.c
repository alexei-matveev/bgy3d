/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fft.h"          /* ComputeFFTfromVec, ... */
#include "hnc3d.h"

typedef struct HNCField
{
  PetscScalar h, c;
} HNCField;

typedef struct HNC3dNewtonStruct
{
  DA da, da1;
  Vec pot;
  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;
  Vec pre;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *c_fft, *h_fft, *ch_fft;

  ProblemData *PD;
} *HNC3dNewtonData;


static void UNUSED_SetBoundaryValue(HNC3dNewtonData HD, Vec g, int x[3], int  n[3])//, real c)
{
  HNCField ***g_vec;
  int i[3], ic[3], k[3], j[3];

  DAVecGetArray(HD->da, g, &g_vec);

  FOR_DIM
    {
      ic[0]=(dim)%3;
      ic[1]=(dim+1)%3;
      ic[2]=(dim+2)%3;
      if( x[ic[0]] == 0 )
	{
	  i[ic[0]] = 0;
	  k[ic[0]] = 1;
	  j[ic[0]] = 2;
	  for(i[ic[1]]= x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++)
	    for(i[ic[2]]=x[ic[2]] ; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++)
	      {
		k[ic[1]]=i[ic[1]];
		k[ic[2]]=i[ic[2]];
		j[ic[1]]=i[ic[1]];
		j[ic[2]]=i[ic[2]];
		g_vec[i[2]][i[1]][i[0]].h = 2*g_vec[k[2]][k[1]][k[0]].h
		  -g_vec[j[2]][j[1]][j[0]].h;

	      }
	}
      if( x[ic[0]]+n[ic[0]] == HD->PD->N[ic[0]])
	{
	  i[ic[0]] = HD->PD->N[ic[0]]-1;
	  k[ic[0]] = HD->PD->N[ic[0]]-2;
	  j[ic[0]] = HD->PD->N[ic[0]]-3;
	  for(i[ic[1]]= x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++)
	    for(i[ic[2]]=x[ic[2]] ; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++)
	      {
		k[ic[1]]=i[ic[1]];
		k[ic[2]]=i[ic[2]];
		j[ic[1]]=i[ic[1]];
		j[ic[2]]=i[ic[2]];
		g_vec[i[2]][i[1]][i[0]].h = 2*g_vec[k[2]][k[1]][k[0]].h
		  -g_vec[j[2]][j[1]][j[0]].h;

	      }
	}
    }

  DAVecRestoreArray(HD->da, g, &g_vec);
}


static HNC3dNewtonData HNC3dNewtonData_malloc(ProblemData *PD)
{
  HNC3dNewtonData HD;
  DA da;
  int n[3], x[3], i[3], N[3];
  PetscScalar  interval[2];
  HNCField ***pot_vec;
  PetscScalar r[3], r_s, L, h[3];
  real **x_M;
  int bufsize, k, N_M;
  real epsilon, sigma;


  HD = (HNC3dNewtonData) malloc(sizeof(*HD));

  HD->LJ_params[0] = 1.0;       /* espilon */
  HD->LJ_params[1] = 1.0;       /* sigma */
  epsilon = HD->LJ_params[0];
  sigma = HD->LJ_params[1];

  HD->beta = PD->beta;
  HD->rho  = PD->rho;


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
	     2,1, PETSC_NULL,PETSC_NULL,PETSC_NULL, &(HD->da));
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	     1,1, PETSC_NULL,PETSC_NULL,PETSC_NULL, &(HD->da1));


  da = HD->da;

  /* Create global vectors */
  DACreateGlobalVector(da, &(HD->pot));
  VecDuplicate(HD->pot, &(HD->pre));

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
  PetscPrintf(PETSC_COMM_WORLD,"HNCNewton ignores solute data from file \"molecule\"!\n");
  if( N_M >1 )
    {
      PetscPrintf(PETSC_COMM_WORLD,"HNCNewton can only handle 1 atom in solute!\n");

    }

  VecSet(HD->pot,0.0);
  DAVecGetArray(da, HD->pot, &pot_vec);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{


	  /* set force vector */
	  /* loop over particles and grid */
	  for(k=0; k<N_M; k++)
	    {
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim]-0.0;
		  if( i[dim] >= N[dim]/2)
		    r[dim]-=L;
		}

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );
	      pot_vec[i[2]][i[1]][i[0]].h +=
		Lennard_Jones( r_s, epsilon, sigma);
	    }

	}



  DAVecRestoreArray(da, HD->pot, &pot_vec);


/*   VecView(HD->pot,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */



  /* Create plan for 3d fft */
  HD->fft_plan = fft_3d_create_plan(PETSC_COMM_WORLD,
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
  if(HD->fft_plan == NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Failed to get fft_plan of proc %d.\n",
		  PD->id);
      exit(1);
    }

  HD->c_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));
  HD->h_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));
  HD->ch_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(FFT_DATA));


  HD->PD=PD;

  Molecule_free(x_M, N_M);

  return HD;
}


static void HNC3dNewtonData_free(HNC3dNewtonData HD)
{
  VecDestroy(HD->pot);
  VecDestroy(HD->pre);
  DADestroy(HD->da);
  DADestroy(HD->da1);
  fft_3d_destroy_plan(HD->fft_plan);

  free(HD->c_fft);
  free(HD->h_fft);
  free(HD->ch_fft);
  free(HD);
}

/* F for Newton method, solving h and c simultaneously */
static PetscErrorCode ComputeHNC_F(SNES snes, Vec g, Vec f, void *pa)
{
  (void) snes;                  /* FIXME: interface obligation? */

  HNC3dNewtonData HD;
  ProblemData *PD;
  DA da;
  int n[3], x[3], i[3], index;
  FFT_DATA *c_fft, *h_fft, *ch_fft;
  HNCField ***g_vec, ***f_vec, ***pot_vec, ***pre_vec;
  PetscScalar c, h, exppot, beta, rho, factor;


  HD = (HNC3dNewtonData) pa;
  da = HD->da;
  PD = HD->PD;
  beta = HD->beta;
  rho = HD->rho;

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  c_fft = HD->c_fft;
  h_fft = HD->h_fft;
  ch_fft = HD->ch_fft;

  DAVecGetArray(da, g, &g_vec);
  DAVecGetArray(da, f, &f_vec);
  DAVecGetArray(da, HD->pot, &pot_vec);
  DAVecGetArray(da, HD->pre, &pre_vec);
  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  c = g_vec[i[2]][i[1]][i[0]].c;
	  h = g_vec[i[2]][i[1]][i[0]].h;
	  exppot = exp(-beta*pot_vec[i[2]][i[1]][i[0]].h+h-c);

	  c_fft[index].re =  c;
	  c_fft[index].im =  0;
	  h_fft[index].re =  h;
	  h_fft[index].im =  0;

	  f_vec[i[2]][i[1]][i[0]].h = h+1-exppot;
	  f_vec[i[2]][i[1]][i[0]].c = c+1-exppot;

	  /* preconditioner */
/* 	  pre_vec[i[2]][i[1]][i[0]].h = 1-exppot; */
/* 	  pre_vec[i[2]][i[1]][i[0]].c = 1+exppot; */

	  index++;
	}

  /* 3d FFT */
  fft_3d(c_fft, c_fft, 1, HD->fft_plan);
  fft_3d(h_fft, h_fft, 1, HD->fft_plan);

  /* convolution */
  for(index=0; index<n[0]*n[1]*n[2]; index++)
    {
      ch_fft[index].re = ( c_fft[index].re*h_fft[index].re
			   -c_fft[index].im*h_fft[index].im );
      ch_fft[index].im = ( c_fft[index].re*h_fft[index].im
			   +c_fft[index].im*h_fft[index].re );
    }

  /* 3d FFT inverse */
  fft_3d(ch_fft, ch_fft, -1, HD->fft_plan);

  index=0;
  factor = rho*PD->h[0]*PD->h[1]*PD->h[2]/PD->N[0]/PD->N[1]/PD->N[2];
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  f_vec[i[2]][i[1]][i[0]].c += factor*ch_fft[index].re;
	  //pre_vec[i[2]][i[1]][i[0]].h = factor*ch_fft[index].re;
	  index++;
	}

  DAVecRestoreArray(da, g, &g_vec);
  DAVecRestoreArray(da, f, &f_vec);
  DAVecRestoreArray(da, HD->pot, &pot_vec);
  DAVecRestoreArray(da, HD->pre, &pre_vec);

/*   VecOutput_hc(HD, g, 0); */
/*   exit(1); */

   //SetBoundaryValue(HD, f, x, n, 0);

  /* preconditioner */
  VecReciprocal(HD->pre);

  return 0;
}

#if PETSC_VERSION_MAJOR > 2
static PetscErrorCode ComputeHNC_Preconditioner (PC pc, Vec x, Vec y)
{
  HNC3dNewtonData HD;
  PCShellGetContext (pc, (void**) &HD);

  PetscErrorCode ierr = VecPointwiseMult (y, HD->pre, x);

  return ierr;
}
#else
static PetscErrorCode ComputeHNC_Preconditioner (void *pa, Vec x, Vec y)
{
  HNC3dNewtonData HD;
  HD = (HNC3dNewtonData) pa;

  PetscErrorCode ierr = VecPointwiseMult (y, HD->pre, x);

  return ierr;
}
#endif

static Vec Compute_gfromhc(HNC3dNewtonData HD, Vec hc)
{
  PetscScalar ***g_vec;
  HNCField ***hc_vec;
  int n[3], x[3], i[3];
  Vec g;

  DACreateGlobalVector(HD->da1, &g);


  DAGetCorners(HD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  DAVecGetArray(HD->da1, g, &g_vec);
  DAVecGetArray(HD->da, hc, &hc_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  g_vec[i[2]][i[1]][i[0]] = hc_vec[i[2]][i[1]][i[0]].h+1;
	}
  DAVecRestoreArray(HD->da1, g, &g_vec);
  DAVecRestoreArray(HD->da, hc, &hc_vec);

  return g;
}

static void UNUSED_VecOutput_hc(HNC3dNewtonData HD, Vec hc, int horc)
{
  PetscScalar ***g_vec;
  HNCField ***hc_vec;
  int n[3], x[3], i[3];
  Vec g;

  DACreateGlobalVector(HD->da1, &g);


  DAGetCorners(HD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  DAVecGetArray(HD->da1, g, &g_vec);
  DAVecGetArray(HD->da, hc, &hc_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  if( horc==0)
	    g_vec[i[2]][i[1]][i[0]] = hc_vec[i[2]][i[1]][i[0]].h;
	  else
	    g_vec[i[2]][i[1]][i[0]] = hc_vec[i[2]][i[1]][i[0]].c;
	}
  DAVecRestoreArray(HD->da1, g, &g_vec);
  DAVecRestoreArray(HD->da, hc, &hc_vec);

  VecView(g,PETSC_VIEWER_STDERR_WORLD);
  VecDestroy(g);

}

static void CreateInitialGuess_HNC(HNC3dNewtonData HD, Vec hc)
{
  HNCField ***hc_vec, ***pot_vec;
  int n[3], x[3], i[3];
  real beta;

  DAGetCorners(HD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  beta = HD->beta;

  DAVecGetArray(HD->da, hc, &hc_vec);
  DAVecGetArray(HD->da, HD->pot, &pot_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  hc_vec[i[2]][i[1]][i[0]].h = exp(-beta*pot_vec[i[2]][i[1]][i[0]].h);
	  hc_vec[i[2]][i[1]][i[0]].c = exp(-beta*pot_vec[i[2]][i[1]][i[0]].h);
	}
  DAVecRestoreArray(HD->da, HD->pot, &pot_vec);
  DAVecRestoreArray(HD->da, hc, &hc_vec);

}

/* solving for h and c of HNC equation with Newton */
static Vec UNUSED_HNC3dNewton_solve(ProblemData *PD, Vec g_ini)
{
  Vec F, hc, g;
  HNC3dNewtonData HD;
  SNES snes;
  KSP ksp;
  PC  pc;
  PetscTruth flg;

  assert(g_ini==PETSC_NULL);
  HD = HNC3dNewtonData_malloc(PD);

  /* Create global vectors */
  VecDuplicate(HD->pot, &F);
  VecDuplicate(HD->pot, &hc);

  /* initial guess */
  //VecSet(hc, 0.0);
  CreateInitialGuess_HNC(HD, hc);

  /* Create the snes environment */
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp, &pc);

  /* set rtol, atol, dtol, maxits */
  // KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType(snes, SNESLS);

  flg = bgy3d_getopt_test ("--user-precond");
  if (flg) { /* user-defined precond */
    /* Set user defined preconditioner */
    PCSetType (pc, PCSHELL);
    PCShellSetContext (pc, HD);
    PCShellSetApply (pc, ComputeHNC_Preconditioner);
  } else
  /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
  PCSetType( pc, PCLU);

/*   ComputeHNC_F(snes, g, F, (void*) HD);  */
/*   VecView(F,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);   */

  SNESSetFunction(snes, F, ComputeHNC_F, HD);

  /* runtime options will override default parameters */
  SNESSetFromOptions(snes);

  /* solve problem */
  SNESSolve(snes, PETSC_NULL, hc);



  /* write out solution */
  SNESGetSolution(snes, &hc);

  g = Compute_gfromhc(HD, hc);


  /* free stuff */
  HNC3dNewtonData_free(HD);

  VecDestroy(F);
  VecDestroy(hc);
  SNESDestroy(snes);

  return g;

}
