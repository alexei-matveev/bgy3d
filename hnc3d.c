/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fft.h"          /* ComputeFFTfromVec, ... */
#include "hnc3d.h"

typedef struct HNC3dDataStruct
{
  DA da;
  Vec pot;
  Vec h_ini;
  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;

  /* things for arbitrary molecule shape */
  Vec c, v;
  FFT_DATA *c_fft, *h_fft, *ch_fft;

  const ProblemData *PD;
} *HNC3dData;

static void Compute_c_HNC(HNC3dData HD, Vec g, Vec c, int x[3], int n[3]);
static PetscErrorCode ComputeHNC2_F(SNES snes, Vec h, Vec f, void *pa);

HNC3dData HNC3dData_malloc(const ProblemData *PD)
{
  HNC3dData HD;
  DA da;
  int n[3], x[3], i[3], N[3];
  PetscScalar ***pot_vec, interval[2], ***hini_vec;
  PetscScalar r[3], r_s, L, h[3], beta;
  real **x_M, h_c1d;
  int bufsize, k, N_M, N_c1d, index;
  PetscViewer pview;
  real epsilon, sigma;


  HD = (HNC3dData) malloc(sizeof(*HD));

  HD->LJ_params[0] = 1.0;       /* epsilon */
  HD->LJ_params[1] = 1.0;       /* sigma */
  epsilon = HD->LJ_params[0];
  sigma = HD->LJ_params[1];

  HD->beta = PD->beta;
  HD->rho  = PD->rho;
  beta = PD->beta;

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
	     1,1, PETSC_NULL,PETSC_NULL,PETSC_NULL, &(HD->da));

  da = HD->da;

  /* Create global vectors */
  DACreateGlobalVector(da, &(HD->pot));
  VecDuplicate(HD->pot, &(HD->c));
  VecDuplicate(HD->pot, &(HD->v));
  VecDuplicate(HD->pot, &(HD->h_ini));

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if( verbosity >2)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Subgrids on processes:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "id %d of %d: %d %d %d\t%d %d %d\n",
			      PD->id, PD->np, x[0], x[1], x[2], n[0], n[1], n[2]);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }


  /*
    Load  c_1d  from  file.   FIXME:  this is  not  needed  for,  say,
    hnc3d_solve()   and  should   probably   be  layed   out  into   a
    subprogram:
  */
  if (1)
    {
      Vec c_1d;
      PetscScalar ***c_vec, *c1d_vec;

      /* c_1d has to be on a grid [0,L] with L=interval[1]-interval[0] */
      PetscViewerBinaryOpen(PETSC_COMM_SELF,"c1dfile",FILE_MODE_READ,
                            &pview);
      VecLoad( pview, VECSEQ, &c_1d);
      VecGetSize(c_1d,&N_c1d);
      h_c1d = L/N_c1d;

      DAVecGetArray(da, HD->c, &c_vec);
      VecGetArray(c_1d, &c1d_vec);

      /* loop over local portion of grid */
      for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
        for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
          for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
            {
              /* set c-vector */
              /* c lives on different grid -> [0,L]^3 (for fft) */
              FOR_DIM
                {
                  r[dim] = i[dim]*h[dim];
                  if( i[dim] >= N[dim]/2 )
                    r[dim] -= L;
                }
              r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );
              index =(int) floor(r_s/h_c1d);
              c_vec [i[2]][i[1]][i[0]] =
                (c1d_vec[index]+ (c1d_vec[index+1]-c1d_vec[index])/h_c1d*
                 (r_s-index*h_c1d));

            }
      DAVecRestoreArray(da, HD->c, &c_vec);
      VecRestoreArray(c_1d, &c1d_vec);
      VecDestroy(c_1d);
    }


  /* Load molecule from file */
  x_M = Load_Molecule(&N_M);



  VecSet(HD->pot,0.0);
  VecSet(HD->h_ini,1.0);
  DAVecGetArray(da, HD->pot, &pot_vec);
  DAVecGetArray(da, HD->h_ini, &hini_vec);

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
		  r[dim] = i[dim]*h[dim]+interval[0]-x_M[k][dim];
/* 		  if( i[dim] >= N[dim]/2 ) */
/* 		    r[dim] -= L; */
		}
	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );
	      pot_vec[i[2]][i[1]][i[0]] +=
		Lennard_Jones( r_s, epsilon, sigma);

	      hini_vec[i[2]][i[1]][i[0]] *=
		exp(-beta* Lennard_Jones( r_s, epsilon, sigma));

	    }
	}

  DAVecRestoreArray(da, HD->h_ini, &hini_vec);
  DAVecRestoreArray(da, HD->pot, &pot_vec);
  VecShift(HD->h_ini, -1.0);

/*    VecView(HD->h_ini,PETSC_VIEWER_STDERR_WORLD);  */
/*    exit(1);   */



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

  HD->c_fft = NULL;
  HD->c_fft = ComputeFFTfromVec(HD->da, HD->fft_plan, HD->c, HD->c_fft);
  HD->h_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2], sizeof(FFT_DATA));
  HD->ch_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2], sizeof(FFT_DATA));


  HD->PD=PD;

  Molecule_free(x_M, N_M);

  return HD;
}


static void HNC3dData_free(HNC3dData HD)
{
  VecDestroy(HD->pot);
  VecDestroy(HD->c);
  VecDestroy(HD->v);
  DADestroy(HD->da);
  fft_3d_destroy_plan(HD->fft_plan);

  free(HD->c_fft);
  free(HD->h_fft);
  free(HD->ch_fft);


  free(HD);
}


static void Compute_cgfft (HNC3dData HD, FFT_DATA *c_fft, FFT_DATA *cg_fft, const int x[3]
		   , const int  n[3], const real h[3])
{
  int i[3], index=0;
  real rho;
  real nenner, re, im, h3;
  Vec t;
  PetscScalar ***t_vec;
  VecDuplicate(HD->pot, &t);

  DAVecGetArray(HD->da, t, &t_vec);
  rho = HD->rho;
  h3 = (h[0]*h[1]*h[2]);

  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* cg_fft = rho*c_fft^2/(1-rho*c_fft) */
	  re = c_fft[index].re*h3;
	  im = c_fft[index].im*h3;
	  nenner = SQR(1.0-rho*re)+SQR(rho*im);

	  cg_fft[index].re = (rho*(SQR(re)-SQR(im))*(1.-rho*re)
			      - 2.*SQR(rho)*re*SQR(im)) / nenner;
	  cg_fft[index].im = (SQR(rho)*(SQR(re)-SQR(im))*im +
			      2.*rho*re*im*(1.-rho*re)) / nenner;

	  //cg_fft[index].im =0;
	  t_vec[i[2]][i[1]][i[0]] = c_fft[index].im;
	  index++;
	}
  DAVecRestoreArray(HD->da, t, &t_vec);
/*   VecView(t,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */
  VecDestroy(t);


}


/* Solve h and c of HNC equation simultaneously, fixpoint iteration */
Vec hnc3d_solve (const ProblemData *PD, Vec g_ini)
{
  HNC3dData HD;
  Vec c, c_old, g, g_old, gg;
  real g_norm, iL3;
  int k, n[3], x[3];
  FFT_DATA *c_fft, *cg_fft;

  assert(g_ini==PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD,"Solving 3d-HNC equation. Fixpoint iteration.\n");

  /* Mixing parameter */
  const real lambda = PD->lambda;

  if(lambda>1 || lambda<0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lambda out of range: lambda=%f\n",lambda);
      exit(1);
    }

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* Convergence threshold: */
  const real norm_tol = PD->norm_tol;

  HD = HNC3dData_malloc(PD);

  DAGetCorners(HD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  iL3 = 1./pow(PD->interval[1]-PD->interval[0],3);

  VecDuplicate(HD->pot, &c);
  VecDuplicate(HD->pot, &g);
  VecDuplicate(HD->pot, &g_old);
  VecDuplicate(HD->pot, &c_old);
  VecDuplicate(HD->pot, &gg);

  /* Set initial guess */
  VecSet(g,0.0);
  VecSet(c,0.0);

  /* set fft data */
  c_fft = NULL;
  cg_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(*cg_fft));

  /* do the iteration */
  for(k=0; k<max_iter; k++)
    {
      /* if(k>3) */
      /*   lambda =0.1; */

      /* set g_old=g */
      VecCopy(g, g_old);
      VecCopy(c, c_old);

      Compute_c_HNC(HD, g, c, x, n);

      /* simple mixing: c = lambda*c+(1-lambda)*c_old */
      VecAXPBY(c, (1-lambda), lambda, c_old);

      c_fft = ComputeFFTfromVec(HD->da, HD->fft_plan, c, c_fft);


      Compute_cgfft(HD, c_fft, cg_fft, x, n, PD->h);

      ComputeVecfromFFT(HD->da, HD->fft_plan, g, cg_fft);

      VecScale(g, iL3);
/*       VecView(c,PETSC_VIEWER_STDERR_WORLD);   */
/*       exit(1);   */
      /* gg=g-g_old */
      VecWAXPY(gg, -1.0, g_old, g);

      VecNorm(gg, NORM_2, &g_norm);
      PetscPrintf(PETSC_COMM_WORLD,"iter %d: norm of difference: %e\t%f\n", k,
		  g_norm, lambda);

      if (g_norm < norm_tol)
        break;
    }

  /* g= gamma+c+1 */
  VecAXPY(g, 1.0, c);
  VecShift(g, 1.0);

  //VecCopy(c,g);

  /* free stuff */
  VecDestroy(c);
  VecDestroy(gg);
  VecDestroy(c_old);
  VecDestroy(g_old);

  free(c_fft);
  free(cg_fft);

  HNC3dData_free(HD);

  bgy3d_save_vec ("g00.bin", g);
  return g;
}


static void Compute_c_HNC(HNC3dData HD, Vec g, Vec c, int x[3], int n[3])
{
  int i[3];
  PetscScalar ***g_vec, ***c_vec, ***pot_vec;
  real beta;

  beta = HD->beta;

  DAVecGetArray(HD->da, g, &g_vec);
  DAVecGetArray(HD->da, c, &c_vec);
  DAVecGetArray(HD->da, HD->pot, &pot_vec);

  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* c=exp(-beta*U+g) */
	  c_vec[i[2]][i[1]][i[0]] = exp(-beta*pot_vec[i[2]][i[1]][i[0]]+
					g_vec[i[2]][i[1]][i[0]]);
	}
  DAVecRestoreArray(HD->da, g, &g_vec);
  DAVecRestoreArray(HD->da, c, &c_vec);
  DAVecRestoreArray(HD->da, HD->pot, &pot_vec);

  /* c=c-g-1 */
  VecAXPY(c, -1.0, g);
  VecShift(c, -1.0);


}


/* solving for h only of HNC equation with Newton */
/* c appears as an input here */
Vec HNC3dNewton2_solve(const ProblemData *PD, Vec g_ini)
{
  Vec F, h;
  HNC3dData HD;
  SNES snes;
  KSP ksp;
  PC  pc;

  assert(g_ini==PETSC_NULL);

  HD = HNC3dData_malloc(PD);

  /* Create global vectors */
  VecDuplicate(HD->pot, &F);
  VecDuplicate(HD->pot, &h);

  /* initial guess */
  VecSet(h, 0.0);
  VecCopy(HD->c, h);

  /* Create the snes environment */
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp, &pc);

  /* set rtol, atol, dtol, maxits */
  // KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType(snes, SNESLS);

  /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
  PCSetType( pc, PCNONE);

/*   ComputeHNC_F(snes, g, F, (void*) HD);  */
/*   VecView(F,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);   */

  SNESSetFunction(snes, F, ComputeHNC2_F, HD);

  /* runtime options will override default parameters */
  SNESSetFromOptions(snes);

  /* solve problem */
  SNESSolve(snes, PETSC_NULL, h);



  /* write out solution */
  SNESGetSolution(snes, &h);




  /* free stuff */
  HNC3dData_free(HD);

  VecDestroy(F);

  SNESDestroy(snes);

  /* g=h+1 */
  VecShift(h,1.0);

  return h;

}

/* F for solving HNC equation with Newton */
/* c appears as an input here */
static PetscErrorCode ComputeHNC2_F(SNES snes, Vec h, Vec f, void *pa)
{
  (void) snes;                  /* FIXME: interface obligation? */

  FFT_DATA *h_fft, *c_fft, *ch_fft;
  int x[3], n[3], i[3], index;
  HNC3dData HD;
  real rho, beta;
  PetscScalar ***f_vec, ***pot_vec, ***v_vec;

  HD = (HNC3dData) pa;
  const ProblemData *PD = HD->PD;

  DAGetCorners(HD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  VecSet(HD->v, 0.0);

  /* fft(h) */
  ComputeFFTfromVec(HD->da, HD->fft_plan, h, HD->h_fft);

  c_fft = HD->c_fft;
  h_fft = HD->h_fft;
  ch_fft= HD->ch_fft;

  rho = HD->rho;
  beta= HD->beta;
  /* fft(h)*fft(c) */
  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  ch_fft[index].re = c_fft[index].re*h_fft[index].re
	                    -c_fft[index].im*h_fft[index].im;
	  ch_fft[index].im = c_fft[index].re*h_fft[index].im
	                    +c_fft[index].im*h_fft[index].re;
	  index++;
	}

  /* v = fft^-1(fft(c)*fft(h)) */
  ComputeVecfromFFT(HD->da, HD->fft_plan, HD->v, ch_fft);

  VecScale(HD->v, PD->h[0]*PD->h[1]*PD->h[2]/PD->N[0]/PD->N[1]/PD->N[2]);





  DAVecGetArray(HD->da, HD->pot, &pot_vec);
  DAVecGetArray(HD->da, HD->v, &v_vec);
  DAVecGetArray(HD->da, f, &f_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  f_vec[i[2]][i[1]][i[0]] = -exp(-beta*pot_vec[i[2]][i[1]][i[0]]+
					 rho*v_vec[i[2]][i[1]][i[0]]);
	}
  DAVecRestoreArray(HD->da, HD->pot, &pot_vec);
  DAVecRestoreArray(HD->da, HD->v, &v_vec);
  DAVecRestoreArray(HD->da, f, &f_vec);

  /* f=f+h+1 */
  VecAXPY(f, 1.0, h);
  VecShift(f, 1.0);

  //VecView(HD->v,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);
  return 0;
}


/* Solving for h of HNC eqauation with fixpoint iteration */
/* c is input */
Vec HNC3d_Solve_h(const ProblemData *PD, Vec g_ini)
{
  HNC3dData HD;
  Vec h, h_old, gg, v; // c
  real g_norm, rho, beta;
  int slow_iter=10, k, n[3], x[3], i[3], index;
  FFT_DATA *c_fft, *ch_fft, *h_fft;
  PetscScalar ***pot_vec, ***h_vec, ***v_vec;

  assert(g_ini==PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD,"Solving 3d-HNC equation. Fixpoint iteration.Fixed c.\n");

  /* Mixing parameter: */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  /* Number of iterations with lambda */
  bgy3d_getopt_int ("--slow-iter", &slow_iter);

  HD = HNC3dData_malloc(PD);
  rho = HD->rho;
  beta = HD->beta;

  DAGetCorners(HD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));




  VecDuplicate(HD->pot, &h);
  VecDuplicate(HD->pot, &h_old);
  VecDuplicate(HD->pot, &gg);
  /* c=HD->c; */
  v=HD->v;
  /* Set initial guess */
  //VecSet(h,0.0);
  VecCopy(HD->h_ini, h);
  VecScale(h, rho);
  //VecShift(h,1);

  /* set fft data */
  c_fft = HD->c_fft;
  ch_fft = HD->ch_fft;
  h_fft = HD->h_fft;

  /* do the iteration */
  for(k=0; k<max_iter; k++)
    {
/*       if(k>slow_iter) */
/* 	lambda =0.9; */



      /* set h_old=h */
      VecCopy(h, h_old);

      /* new */
      //VecShift(h,-1);
      /* fft(h) */
      ComputeFFTfromVec(HD->da, HD->fft_plan, h, h_fft);

      /* fft(h)*fft(c) */
      index=0;
      /* set int(h)=0 for numerical stabilization */
      h_fft[0].re=0;h_fft[0].im=0;
      /* loop over local portion of grid */
      for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
	for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
	  for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	    {
	      ch_fft[index].re = c_fft[index].re*h_fft[index].re
		-c_fft[index].im*h_fft[index].im;
	      ch_fft[index].im = c_fft[index].re*h_fft[index].im
		+c_fft[index].im*h_fft[index].re;
	      index++;
	    }

      /* v=fft^-1(fft(h)*fft(c)) */
      ComputeVecfromFFT(HD->da, HD->fft_plan, v, ch_fft);

      VecScale(v, PD->h[0]*PD->h[1]*PD->h[2]/PD->N[0]/PD->N[1]/PD->N[2]);

      DAVecGetArray(HD->da, HD->pot, &pot_vec);
      DAVecGetArray(HD->da, v, &v_vec);
      DAVecGetArray(HD->da, h, &h_vec);

      /* loop over local portion of grid */
      for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
	for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
	  for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	    {
	      h_vec[i[2]][i[1]][i[0]]=rho*(exp(-beta*pot_vec[i[2]][i[1]][i[0]]
					       + v_vec[i[2]][i[1]][i[0]])-1.0);
/* 	      h_vec[i[2]][i[1]][i[0]]=(exp(-beta*pot_vec[i[2]][i[1]][i[0]] */
/* 					   + rho*v_vec[i[2]][i[1]][i[0]])-0.0); */
	    }
/*       PetscPrintf(PETSC_COMM_WORLD,"%e\t%e\t%e\n", h_vec[0][0][0], pot_vec[0][0][0],  */
/*       		  v_vec[0][0][0]); */

      DAVecRestoreArray(HD->da, HD->pot, &pot_vec);

/*       DAVecGetArray(HD->da, h_old, &pot_vec); */
/*       PetscPrintf(PETSC_COMM_WORLD,"%e\t%e\n",h_vec[0][64][64], h_vec[0][64][64]-pot_vec[0][64][64]); */
/*       g_norm = h_vec[0][64][64]-pot_vec[0][64][64]; */
/*       DAVecRestoreArray(HD->da, h_old, &pot_vec); */
/*       VecShift(h,-g_norm); */

      DAVecRestoreArray(HD->da, v, &v_vec);
      DAVecRestoreArray(HD->da, h, &h_vec);



      /* gg=h-h_old */
      VecWAXPY(gg, -1.0, h_old, h);
      VecNorm(gg, NORM_INFINITY, &g_norm);
      PetscPrintf(PETSC_COMM_WORLD,"iter %d: norm of difference: %e\t%f\n",
		  k+1, g_norm, lambda);

      if(g_norm < norm_tol)
	{
	  /* simple mixing: h = lambda*h+(1-lambda)*h_old */
	  VecAXPBY(h, (1-lambda), lambda, h_old);
	  break;
	}
      else
	/* simple mixing: h = lambda*h+(1-lambda)*h_old */
	VecAXPBY(h, (1-lambda), lambda, h_old);

/*       VecSum(h, &g_norm); */
/*       VecSum(h_old, &gold_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD,"%e\n",g_norm*PD->h[0]*PD->h[1]*PD->h[2]/1000); */
/*       VecShift(h,-g_norm*PD->h[0]*PD->h[1]*PD->h[2]/1000); */

    }

  /* g= h+1 */
  VecScale(h,1./rho);
  VecShift(h, 1.0);

  //VecCopy(c,g);

  /* free stuff */
  VecDestroy(h_old);
  VecDestroy(gg);


  HNC3dData_free(HD);

  return h;
}
