/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2OS.c,v 1.20 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"

#include "bgy3d-solvents.h"
#include "bgy3d-getopt.h"
#include "bgy3d-solutes.h"
#include "bgy3dH2O.h"
#include "bgy3dH2OS.h"
#include "bgy3d-fft.h"

#ifndef L_BOUNDARY_MG
#include "bgy3d-multigrid.h"    /* InitializeDMMGSolver */
#endif

#ifdef WITH_COMPLEX
#include <complex.h>
#endif

#define damp0 0.001

extern real NORM_REG;

/*
 * These are  the two  solvent sites.  Coordinates  will not  be used.
 * Respective parameters are #defined  elsewhere. Also do not take the
 * names  of the  sites  literally.   The same  structure  is used  to
 * represent all (2-site) solvents, such as HCl.
 *
 * FIXME:  at many places  it is  assumed that  the number  of solvent
 * sites is exactly two:
 */
static const Site solvent[] =
  {{"h", {0.0, 0.0, 0.0}, sH, eH, qH}, /* dont use sH, eH, qH below */
   {"o", {0.0, 0.0, 0.0}, sO, eO, qO}}; /* same for sO, eO, qO */

static void load (const State *BHD, Vec g2[2][2]);

static State initialize_state (const ProblemData *PD)
{
  State BHD;
  PetscErrorCode ierr;

  /****************************************************/
  /* set Lennard-Jones and Coulomb parameters */
  /****************************************************/

  /* water hydrogen */
  BHD.LJ_paramsH[0] = eH;  /* epsilon  */
  BHD.LJ_paramsH[1] = sH;  /* sigma    */
  BHD.LJ_paramsH[2] = SQR(qH); /* charge product */

  /* water oxygen */
  BHD.LJ_paramsO[0] = eO;  /* epsilon  */
  BHD.LJ_paramsO[1] = sO;  /* sigma    */
  BHD.LJ_paramsO[2] = SQR(qO); /* charge product */

  /* water O-H mixed parameters */
  BHD.LJ_paramsHO[0] = sqrt(eH*eO);  /* epsilon  */
  BHD.LJ_paramsHO[1] = 0.5*(sH+sO);  /* sigma    */
  BHD.LJ_paramsHO[2] = qH*qO; /* charge product */

  /****************************************************/

  BHD.PD = PD;

  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  //PetscPrintf(PETSC_COMM_WORLD, "Boundary smoothing parameters : SL= %f  SR= %f\n", SL, SR);
  //PetscPrintf(PETSC_COMM_WORLD, "ZEROPAD= %f\n", ZEROPAD);
  //PetscPrintf(PETSC_COMM_WORLD, "Regularization of normalization: NORM_REG= %e\n", NORM_REG);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", PD->beta);
  /******************************/
  BHD.beta = PD->beta;
  BHD.rho  = PD->rho;
  BHD.rhos[0] = PD->rho;
  BHD.rhos[1] = PD->rho;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
#ifndef L_BOUNDARY_MG
  bgy3d_fft_init_da (PD->N, &BHD.fft_plan_fw, &BHD.fft_plan_bw, &BHD.da, NULL);
#else
  /* multigrid, apparently needs two descriptors: */
  bgy3d_fft_init_da (PD->N, &BHD.fft_plan_fw, &BHD.fft_plan_bw, &BHD.da, &BHD.da_dmmg);
#endif

  const DA da = BHD.da;         /* shorter alias */

  /* Create global vectors */
  ierr = DACreateGlobalVector (da, &BHD.g_ini[0]); assert (!ierr);
  ierr = DACreateGlobalVector (da, &BHD.g_ini[1]); assert (!ierr);
  ierr = DACreateGlobalVector (da, &BHD.gHO_ini); assert (!ierr);
  ierr = DACreateGlobalVector (da, &BHD.u2[0][0]); assert (!ierr);
  ierr = DACreateGlobalVector (da, &BHD.u2[1][1]); assert (!ierr);
  ierr = DACreateGlobalVector (da, &BHD.u2[0][1]); assert (!ierr);

  /* Pair quantities  here, use symmetry wrt  (i <-> j)  to save space
     and work: */
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j <= i; j++) {
      /* FIXME: see  bgy3d_load_vec below, arent  we overwriting these
         g2 vecs? */
      ierr = DACreateGlobalVector (da, &BHD.g2[i][j]);
      BHD.g2[j][i] = BHD.g2[i][j];
      assert (!ierr);

      FOR_DIM
        {
          ierr = DACreateGlobalVector (da, &BHD.F[i][j][dim]);
          BHD.F[j][i][dim] = BHD.F[i][j][dim];
          assert (!ierr);

          ierr = DACreateGlobalVector (da, &BHD.F_l[i][j][dim]);
          BHD.F_l[j][i][dim] = BHD.F_l[i][j][dim];
          assert (!ierr);
        }
    }
  }

  ierr = DACreateGlobalVector (da, &BHD.pre); assert (!ierr);

  FOR_DIM
    {
      ierr = DACreateGlobalVector (da, &BHD.v[dim]); assert (!ierr);
    }


  VecSet(BHD.pre,0.0);

#ifdef L_BOUNDARY
   /* Create Matrix with appropriate non-zero structure */
  DAGetMatrix( da, MATMPIAIJ, &BHD.M);
  DACreateGlobalVector(da, &BHD.x_lapl[0]);
  DACreateGlobalVector(da, &BHD.x_lapl[1]);
  VecSet(BHD.x_lapl[0], 0.0);
  VecSet(BHD.x_lapl[1], 0.0);
#endif

  /* Allocate memory for fft */
  FOR_DIM
    {
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j <= i; j++) {
          BHD.f_g2_fft[i][j][dim] = bgy3d_fft_malloc (da);
          BHD.f_g2_fft[j][i][dim] = BHD.f_g2_fft[i][j][dim];

          BHD.fl_g2_fft[i][j][dim] = bgy3d_fft_malloc (da);
          BHD.fl_g2_fft[j][i][dim] = BHD.fl_g2_fft[i][j][dim];
        }
      }

      BHD.fg2_fft[dim] = bgy3d_fft_malloc (da);

      BHD.fO_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fH_fft[dim] = bgy3d_fft_malloc (da);
    }

  BHD.g_fft = bgy3d_fft_malloc (da);
  BHD.gfg2_fft = bgy3d_fft_malloc (da);
  BHD.fft_scratch = bgy3d_fft_malloc (da);
  BHD.u2_fft[0][0] = bgy3d_fft_malloc (da);
  BHD.u2_fft[1][1] = bgy3d_fft_malloc (da);
  BHD.u2_fft[0][1] = bgy3d_fft_malloc (da);
  BHD.wHO_fft = bgy3d_fft_malloc (da);
  BHD.wHH_fft = bgy3d_fft_malloc (da);

  /* FIXME: broken, see the comments inside the function itself: */
  load (&BHD, BHD.g2);

  return BHD;
}


/* This returns  an array of newly allocated  Vecs with solven-solvent
   pair distributions:  */
static void load (const State *BHD, Vec g2[2][2])
{
  /* FIXME: this is broken. In one case (ifdef CS2) it is assumed that
     g2 vectors are  allocated and here they are  filled with data. In
     another case the vectors are constructed from scratch, eventually
     overwriting  the references  to other  allocated Vecs  and making
     them inaccesible. */
#ifdef CS2
  ReadPairDistribution (BHD, "g2C", g2[1][1]);
  ReadPairDistribution (BHD, "g2S", g2[0][0]);
  ReadPairDistribution (BHD, "g2CS", g2[0][1]);
#else
  /* Read g^2 from file. FIXME: why did we allocate them above? */
  g2[0][0] = bgy3d_load_vec ("g00.bin");
  g2[0][1] = bgy3d_load_vec ("g01.bin");
  g2[1][1] = bgy3d_load_vec ("g11.bin");
#endif
  g2[1][0] = g2[0][1];
}


static void finalize_state (State *BHD)
{
  MPI_Barrier( PETSC_COMM_WORLD);

  /* Pair quantities here: */
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j <= i; j++) {
      VecDestroy (BHD->g2[i][j]);
      FOR_DIM
        {
          VecDestroy (BHD->F[i][j][dim]);
          VecDestroy (BHD->F_l[i][j][dim]);

          bgy3d_fft_free (BHD->f_g2_fft[i][j][dim]);
          bgy3d_fft_free (BHD->fl_g2_fft[i][j][dim]);
        }
    }
  }

  FOR_DIM
    {
      VecDestroy(BHD->v[dim]);
      bgy3d_fft_free (BHD->fg2_fft[dim]);

      bgy3d_fft_free (BHD->fO_fft[dim]);
      bgy3d_fft_free (BHD->fH_fft[dim]);
    }
  bgy3d_fft_free (BHD->g_fft);
  bgy3d_fft_free (BHD->gfg2_fft);
  bgy3d_fft_free (BHD->fft_scratch);
  bgy3d_fft_free (BHD->u2_fft[0][0]);
  bgy3d_fft_free (BHD->u2_fft[1][1]);
  bgy3d_fft_free (BHD->u2_fft[0][1]);
  bgy3d_fft_free (BHD->wHO_fft);
  bgy3d_fft_free (BHD->wHH_fft);

  VecDestroy(BHD->g_ini[0]);
  VecDestroy(BHD->g_ini[1]);
  VecDestroy(BHD->gHO_ini);
  VecDestroy(BHD->u2[0][0]);
  VecDestroy(BHD->u2[1][1]);
  VecDestroy(BHD->u2[0][1]);
  VecDestroy(BHD->pre);
#ifdef L_BOUNDARY
  MatDestroy(BHD->M);
  KSPDestroy(BHD->ksp);
  VecDestroy(BHD->x_lapl[1]);
  VecDestroy(BHD->x_lapl[0]);
#endif
#ifdef L_BOUNDARY_MG
  DMMGDestroy(BHD->dmmg);
#endif
  DADestroy(BHD->da);

  fftwnd_mpi_destroy_plan(BHD->fft_plan_fw);
  fftwnd_mpi_destroy_plan(BHD->fft_plan_bw);
}

#ifdef L_BOUNDARY
/* Initialize M-Matrix with appropriate stencil */
void InitializeLaplaceMatrix(State *BHD, real zpad)
{
  Mat M;
  DA da;
  int x[3], n[3], i[3], N[3], border;
  MatStencil col[3],row;
  PetscScalar v[3], vb=1.0;
  real h[3];

  PetscPrintf(PETSC_COMM_WORLD,"Assembling Matrix...");

  da = BHD->da;
  const ProblemData *PD = BHD->PD;
  M = BHD->M;
  MatZeroEntries(M);


  border = (int) ceil( ((PD->interval[1]-PD->interval[0])-(2.*zpad))/PD->h[0]/2. );


  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          FOR_DIM
            {
              col[dim].i=i[0];
              col[dim].j=i[1];
              col[dim].k=i[2];
              row.i=i[0];
              row.j=i[1];
              row.k=i[2];
            }

          /* Boundary */
          if( i[0] <= border+1 || i[1] <= border+1 || i[2] <= border+1)
            MatSetValuesStencil(M,1,&row,1,col+1,&vb,ADD_VALUES);
          else if( i[0]>=N[0]-1-border || i[1]>=N[1]-1-border || i[2]>=N[2]-1-border)
            MatSetValuesStencil(M,1,&row,1,col+1,&vb,ADD_VALUES);
          else
            {
              FOR_DIM
                {
                  /* position in matrix */
                  switch(dim)
                    {
                    case 0: col[0].i -= 1; col[2].i += 1;break;
                    case 1: col[0].j -= 1; col[2].j += 1;break;
                    case 2: col[0].k -= 1; col[2].k += 1;break;
                    }
                  /* values to enter */
                  v[0]=1.0/SQR(h[dim]);
                  v[1]=-2.0/SQR(h[dim]);
                  v[2]=+1.0/SQR(h[dim]);


                  MatSetValuesStencil(M,1,&row,3,col,v,ADD_VALUES);
                  switch(dim)
                    {
                    case 0: col[0].i += 1; col[2].i -= 1;break;
                    case 1: col[0].j += 1; col[2].j -= 1;break;
                    case 2: col[0].k += 1; col[2].k -= 1;break;
                    }

                }
            }
        }


  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);

  PetscPrintf(PETSC_COMM_WORLD,"done.\n");
}

void InitializeKSPSolver(State *BHD)
{
  PC pc;



  /* Create ksp environment */
  KSPCreate( PETSC_COMM_WORLD, &BHD->ksp);
  KSPGetPC(BHD->ksp, &pc);
  KSPSetTolerances(BHD->ksp, 1.0e-4, 1.0e-4, 1.0e+5, 1000);

  /* Set Matrix */
  //KSPSetOperators( BHD->ksp, BHD->M, BHD->M, SAME_NONZERO_PATTERN);
  KSPSetOperators( BHD->ksp, BHD->M, BHD->M, SAME_PRECONDITIONER);
  /* Set preconditioner */
  PCSetType( pc, PCBJACOBI);

  KSPSetInitialGuessNonzero(BHD->ksp, PETSC_TRUE);

  /* runtime options will override default parameters */
  //KSPSetFromOptions(BHD->ksp);

}

static void CopyBoundary (const State *BHD, Vec gfrom, Vec gto, real zpad)
{
  DA da;
  int x[3], n[3], i[3], N[3], border; // ic[3], k;
  PetscScalar ***gfrom_vec, ***gto_vec;

  const ProblemData *PD = BHD->PD;
  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  VecSet(gto, 0.0);

  border = (int) ceil( ((PD->interval[1]-PD->interval[0])-(2.*zpad))/PD->h[0]/2. );

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  DAVecGetArray(da, gfrom, &gfrom_vec);
  DAVecGetArray(da, gto, &gto_vec);


/*   FOR_DIM */
/*     { */
/*       ic[0]=dim; */
/*       ic[1]=(dim+1)%3; */
/*       ic[2]=(dim+2)%3; */
/*       i[ic[0]] = 0; */
/*       if( x[ic[0]] == 0 )  */
/*      for( i[ic[1]]=x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++) */
/*        for( i[ic[2]]=x[ic[2]]; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++) */
/*          gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]]; */
/*       i[ic[0]] = N[ic[0]]-1; */
/*       if( x[ic[0]]+n[ic[0]] == N[ic[0]] )  */
/*      for( i[ic[1]]=x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++) */
/*        for( i[ic[2]]=x[ic[2]]; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++) */
/*          gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]]; */
/*     } */

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          if( ( i[0]<=border+1 || i[0]>=N[0]-1-border) ||
              ( i[1]<=border+1 || i[1]>=N[1]-1-border) ||
              ( i[2]<=border+1 || i[2]>=N[2]-1-border) )
            gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]];
        }


  DAVecRestoreArray(da, gfrom, &gfrom_vec);
  DAVecRestoreArray(da, gto, &gto_vec);
}

real ImposeLaplaceBoundary (const State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter)
{
  real mpi_start, mpi_stop;
  static int count=0;

  count++;

  /* computation time measurement start point */
  MPI_Barrier( PETSC_COMM_WORLD);
  mpi_start = MPI_Wtime();

  /* Get boundary of g */
  CopyBoundary(BHD, g, b, zpad);
  //VecSet(x, 0.0);


  /* Solve Laplace */
  KSPSolve(BHD->ksp, b, x);

  if(iter!=NULL)
    KSPGetIterationNumber(BHD->ksp, iter);


  /* subtract solution from g */
  VecAXPY(g, -1.0, x);

  /* computation time measurement stop point */
  MPI_Barrier( PETSC_COMM_WORLD);
  mpi_stop = MPI_Wtime();

  return mpi_stop-mpi_start;

  /*  VecView(b,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);  */

}


#endif

void ReadPairDistribution (const State *BHD, const char *filename, Vec g2)
{
  DA da;
  FILE *fp;
  real *xg, *g;
  real r[3], r_s, h[3], interval[2];
  int index=0;
  int x[3], n[3], i[3], k;
  PetscScalar ***g2_vec;

  da = BHD->da;
  const ProblemData *PD = BHD->PD;
  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  /* read file */
  fp = fopen(filename, "r");
  if(fp==NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Could not open file %s.\n", filename);
      exit(1);
    }
  xg= (real*) malloc(sizeof(*xg));
  g= (real*) malloc(sizeof(*g));

  while( fscanf(fp,"%lf %lf", &xg[index], &g[index]) == 2)
    {
      index++;
      xg= (real*) realloc(xg, (index+1)*sizeof(*xg));
      g= (real*) realloc(g, (index+1)*sizeof(*g));
    }
  index--;
  PetscPrintf(PETSC_COMM_WORLD,"Read %d lines from file %s, x=[%f,%f]\n", index+1, filename,
              xg[0], xg[index]);
  fclose(fp);


  /* interpolate to 3d grid */

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  DAVecGetArray(da, g2, &g2_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          FOR_DIM
            r[dim] = i[dim]*h[dim]+interval[0];
          r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

          /* find x in array */
          for(k=0; k<=index; k++)
            if(r_s<xg[k])
              break;
          if(k==0)
            g2_vec[i[2]][i[1]][i[0]] = 0;
          else if(k>=index)
            g2_vec[i[2]][i[1]][i[0]] = 1.0;
          else
            g2_vec[i[2]][i[1]][i[0]] = g[k] + (r_s-xg[k])*(g[k+1]-g[k])/(xg[k+1]-xg[k]);
          if(g2_vec[i[2]][i[1]][i[0]]<0)
            g2_vec[i[2]][i[1]][i[0]]=0;

        }
  DAVecRestoreArray(da, g2, &g2_vec);

  free(xg);
  free(g);

}

static void  pair (State *BHD,
                   const real LJ_params[3],
                   const Vec g2,
                   Vec F[3], Vec F_l[3],         /* intent(out) */
                   fftw_complex *f_g2_fft[3], /* intent(out) */
                   fftw_complex *fl_g2_fft[3], /* intent(out) */
                   Vec u2, fftw_complex *u2_fft, /* intent(out) */
                   real damp, real damp_LJ)
{
  PetscScalar ***f_vec[3];
  PetscScalar ***fl_vec[3];
  real r[3], r_s, h[3], interval[2];
  int x[3], n[3], i[3];

  const real epsilon = LJ_params[0]; /* geometric average of two */
  const real sigma = LJ_params[1];   /* arithmetic average of two */
  const real q2 = LJ_params[2];      /* charge product */

  const ProblemData *PD = BHD->PD;
  const DA da = BHD->da;

  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  FOR_DIM
    {
      VecSet (F[dim], 0.0);
      VecSet (F_l[dim], 0.0);
    }

  /* Compute Coulomb from fft part */
  /*   ComputeFFTfromCoulombII(BHD, F, F_l, u2_fft, LJ_params, damp); */

  /* Here u2 is intent(out) in the next call: */
  ComputeFFTfromCoulomb(BHD, u2, F_l, u2_fft, q2, damp0);

  FOR_DIM
    {
      DAVecGetArray (da, F[dim], &f_vec[dim]);
      DAVecGetArray (da, F_l[dim], &fl_vec[dim]);
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

          FOR_DIM
            {
              /* Lennard-Jones */
              f_vec[dim][i[2]][i[1]][i[0]] +=
                damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilon, sigma);

              /* Coulomb short */
              f_vec[dim][i[2]][i[1]][i[0]] +=
                damp * Coulomb_short_grad( r_s, r[dim], q2);
            }
        }

  FOR_DIM
    {
      DAVecRestoreArray (da, F[dim], &f_vec[dim]);
      DAVecRestoreArray (da, F_l[dim], &fl_vec[dim]);
    }

  /* Compute FFT(F * g^2).

     F * g2 = (F_LJ + F_coulomb_short) * g2
            + (F_coulomb_long * g2 - F_coulomb_long)
            + F_coulomb_long

     see  (5.101) and (5.102)  in Jager's  thesis. FFT(F_coulomb_long)
     has been  calculated as u2_fft  by ComputeFFTfromCoulomb() above.
     The code needs at least one work vector, use this: */
  Vec work = BHD->v[0];

  FOR_DIM
    {
      PetscErrorCode err;

      /* First (F_LJ + F_coulomb_short) * g2: */
      err = VecPointwiseMult (work, g2, F[dim]);
      assert (!err);

      /* Next FFT((F_LJ + F_coulomb_short) * g2): */
      ComputeFFTfromVec_fftw (da, BHD->fft_plan_fw,
                              work, f_g2_fft[dim],
                              BHD->fft_scratch);

      /* Now Coulomb long. F_coulomb_long * g2: */
      err = VecPointwiseMult (work, g2, F_l[dim]);
      assert (!err);

      /* Next F_coulomb_long * g2 - F_coulomb_long: */
      VecAXPY(work, -1.0, F_l[dim]);

      /* Finally FFT(F_coulomb_long * g2 - F_coulomb_long): */
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw,
                             work, fl_g2_fft[dim],
                             BHD->fft_scratch);
    }
}

void RecomputeInitialFFTs (State *BHD, real damp, real damp_LJ)
{
  real ff_params[3];

  PetscPrintf (PETSC_COMM_WORLD,
               "Recomputing FFT data with damping factor %f (damp_LJ=%f)\n",
               damp, damp_LJ);

  /* Over all distinct solvent site pairs: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        /* For   debug  only,   check   the  symmetry   of  the   pair
           quantities: */
        if (1)
          {
            assert (BHD->g2[j][i] == BHD->g2[i][j]);
            FOR_DIM
              {
                assert (BHD->F[j][i][dim] == BHD->F[i][j][dim]);
                assert (BHD->F_l[j][i][dim] == BHD->F_l[i][j][dim]);
                assert (BHD->f_g2_fft[j][i][dim] == BHD->f_g2_fft[i][j][dim]);
                assert (BHD->fl_g2_fft[j][i][dim] == BHD->fl_g2_fft[i][j][dim]);
              }
          }

        /* Pair  interaction parameters.   FIXME: ignores  or, rather,
           recomputes BHD->LJ_params* here: */
        ff_params[0] = sqrt (solvent[i].epsilon * solvent[j].epsilon);
        ff_params[1] = 0.5 * (solvent[i].sigma + solvent[j].sigma);
        ff_params[2] = solvent[i].charge * solvent[j].charge;

        /* Does real work: */
        pair (BHD, ff_params,
              BHD->g2[i][j], BHD->F[i][j], BHD->F_l[i][j],
              BHD->f_g2_fft[i][j], BHD->fl_g2_fft[i][j],
              BHD->u2[j][i], BHD->u2_fft[j][i], /* FIXME: ij ji */
              damp, damp_LJ);
      }
}

/*
 * This is supposed to compute the k-grid representation of the BGY3dM
 * equation:
 *
 *     Δu = - β ρ div((F g2) * g)
 *
 * with "*" being  the convolution that in Fourier  rep degenerates to
 * pointwise  product,  g2 and  g  being  the (fixed)  solvent-solvent
 * site-site  distribution  function and  (the  unknown) solvent  site
 * distribution around the solute.
 *
 * The kernel dfg(kx,  ky, kz) includes the Poisson factor  1 / k2 but
 * does not include the facotor (rho * beta):
 *
 *     u(kx, ky, kz) = β ρ dfg(kx, ky, kz) g(kx, ky, kz)
 *
 * The three fg arrays are intent(in). The dfg array is intent(out).
 *
 * No known side effects.
 */
static void kernel (const DA da, const ProblemData *PD, fftw_complex *(fg[3]),
                    const fftw_complex *coul,
                    fftw_complex dfg[])
{
  int x[3], n[3], i[3], N[3], ic[3];

  FOR_DIM
    N[dim] = PD->N[dim];

  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];
  const real L3 = L * L * L;
  const real fac = L / (2.0 * M_PI); /* BHD->f ist nur grad U, nicht F=-grad U  */
  const real scale = fac / L3;

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* Loop over local portion of grid: */
  int ijk = 0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          FOR_DIM
            {
              if( i[dim] <= N[dim]/2)
                ic[dim] = i[dim];
              else
                ic[dim] = i[dim] - N[dim];
            }

          /* phase shift factor for x=x+L/2 */
          real sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);

          real k2 = SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0]);

          real k_fac;
          if (k2 > 0)
              k_fac = sign * h3 * h3 * scale / k2;
          else
              k_fac = 0.0;

          /*
           * Compute the  (Fourier transform of) of  the divergence of
           * the "weighted force"  vector div (F g) which  serves as a
           * convolution  kernel  in BGY3dM  equations.  Note the  the
           * derivative in momentum space involves a factor -i so that
           * real- and imaginary  parts of divergence are proportional
           * to
           *
           *     Re = + dot(k, Im(F g))
           *     Im = - dot(k, Re(F g))
           *
           * The  additional factor  -k^(-2)  effectively included  in
           * k_fac recovers the Fourier transform of the corresponding
           * Poisson solution.
           */

          real re = 0.0;
          real im = 0.0;
          for (int p = 0; p < 3; p++) {
              re += ic[p] * fg[p][ijk].im;
              im -= ic[p] * fg[p][ijk].re;
          }
          dfg[ijk].re = k_fac * re;
          dfg[ijk].im = k_fac * im;

          /*
           * FIXME: Origin of  this occasional addition needs some
           * explanation.
           *
           * The   difference   between   Compute_H2O_interS_C()   and
           * Compute_H2O_interS() of the original code reduces to this
           * addendum.   Supply   a   NULL   pointer   for   coul   in
           * Compute_H2O_interS_C()   to    get   the   behaviour   of
           * Compute_H2O_interS().
           *
           * Long range Coulomb part (right one):
           */
          if (coul) {
              dfg[ijk].re += (h3 / L3) * sign * coul[ijk].re;
              dfg[ijk].im += (h3 / L3) * sign * coul[ijk].im;
          }
          ijk++;
        }
}

#ifdef WITH_COMPLEX
static void mul (int n, const complex *restrict a, const complex *restrict b,
                 const double alpha, complex *restrict c)
{
  for (int i = 0; i < n; i++)
      c[i] += alpha * a[i] * b[i];
}
#endif

/*
 * This applies the kernel compured by  kernel() to FFT of g to obtain
 * an increment to "dg". The latter probably needs a better name. Dont
 * forget to clear "dg" early enough.
 *
 * Complex array dg is intent(inout).
 */
static void apply (const DA da,
                   const fftw_complex *restrict ker,  /* kernel */
                   const fftw_complex *restrict g,    /* current g */
                   const real scale,          /* overall scale */
                   fftw_complex *restrict dg) /* intent(inout) */
{
  int x[3], n[3];

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  int n3 = n[0] * n[1] * n[2];

#ifdef WITH_COMPLEX
  mul (n3, (complex*) ker, (complex*) g, scale, (complex*) dg);
#else
  /* loop over local portion of grid */
  for (int ijk = 0; ijk < n3; ijk++) {
      /*
       * Retrive  the precomuted Fourier  transform of  of the
       * divergence of  the "weighted force" vector  div (F g)
       * which  serves  as  a  convolution  kernel  in  BGY3dM
       * equations.  The additional factor -k^(-2) effectively
       * included in the kernel recovers the Fourier transform
       * of  the corresponding  Poisson solution.   The factor
       * scale =  βρ is not  included, on the other  hand. See
       * kernel() for details:
       */
      dg[ijk].re += scale * (ker[ijk].re * g[ijk].re -
                             ker[ijk].im * g[ijk].im);
      dg[ijk].im += scale * (ker[ijk].re * g[ijk].im +
                             ker[ijk].im * g[ijk].re);
  }
#endif
}

/*
 * Side effects:
 *
 *     Uses BHD->{g_fft, gfg2_fft, fft_scratch} as work arrays.
 */
static void Compute_H2O_interS_C (const State *BHD,
                                  fftw_complex *(fg2_fft[3]), Vec g,
                                  const fftw_complex *coul_fft,
                                  real rho, Vec dg_help)
{
  /* Avoid separate VecScale at the end: */
  const real scale = rho * BHD->PD->beta;

  /************************************************/
  /* rho*F*g^2 g*/
  /************************************************/


  /* fft(g) */
  ComputeFFTfromVec_fftw (BHD->da, BHD->fft_plan_fw, g,
                          BHD->g_fft,        /* result */
                          BHD->fft_scratch); /* work array */

  /* FIXME:  Move  computation  of   the  kernel  out  of  the  BGY3dM
     iterations.   Here  we   put   the  fft   of   the  kernel   into
     BHD->fft_scratch: */
  kernel (BHD->da, BHD->PD, fg2_fft, coul_fft,
          BHD->fft_scratch);    /* result */

  /* Apply the kernel, Put result into BHD->gfg2_fft */
  bgy3d_fft_set (BHD->da, BHD->gfg2_fft, 0.0);

  apply (BHD->da, BHD->fft_scratch, BHD->g_fft, scale,
         BHD->gfg2_fft);        /* will be incremented */

  /* ifft(dg) */
  ComputeVecfromFFT_fftw (BHD->da, BHD->fft_plan_bw,
                          dg_help, /* result */
                          BHD->gfg2_fft,
                          BHD->fft_scratch);
}

/*
 * Side effects:
 *
 *     Uses BHD->{g_fft, gfg2_fft, fft_scratch} as work arrays.
 */
void Compute_H2O_interS (const State *BHD, /* NOTE: modifies BHD->fft dynamic arrays */
                         fftw_complex *(fg2_fft[3]), Vec g,
                         real rho, Vec dg_help)
{
    Compute_H2O_interS_C(BHD, fg2_fft, g, NULL, rho, dg_help);
}

real ComputeCharge(State *BHD, Vec g1, Vec g2)
{
  real g1_sum, g2_sum, c;
  Vec help;


  const ProblemData *PD = BHD->PD;
  help = BHD->v[0];

  VecCopy(g1, help);
  VecShift(help, -1.0);
  VecSum(help, &g1_sum);
  VecCopy(g2, help);
  VecShift(help, -1.0);
  VecSum(g2, &g2_sum);

  c= PD->h[0]*PD->h[1]*PD->h[2]*BHD->rho*(g1_sum-0.*g2_sum);

  return c;
}

typedef struct StepData
{
  State *BHD ;
  Vec dg_newO, dgO, dgH;
} StepData;

static PetscErrorCode ComputeStepFunction(SNES snes, Vec x, Vec f, void *data)
{
  (void) snes;                 /* unused, but interface obligation? */

  StepData *SD;
  real con, sumO, sumH; // res
  PetscScalar *x_vec, *f_vec;
  State *BHD;
  Vec gO, gH, dgO2;

  SD = (StepData*) data;
  BHD = SD->BHD;
  const ProblemData *PD = BHD->PD;
  gO = BHD->v[0];
  gH = BHD->v[1];
  dgO2 = BHD->v[3];


  con = PD->h[0]*PD->h[1]*PD->h[2]*BHD->rho;

  VecCopy(SD->dgO, dgO2);
  VecGetArray(x, &x_vec);
  //fprintf(stderr,"x = %e\n", x_vec[0]);
  // VecAXPY(dgO2, x_vec[0], SD->dg_newO);
  VecAXPBY(dgO2, SQR(x_vec[0]), (1-SQR(x_vec[0])), SD->dg_newO);
  //VecScale(dgO2, x_vec[0]);
  VecRestoreArray(x, &x_vec);

  ComputeH2O_g( gH, BHD->g_ini[0], SD->dgH);
  ComputeH2O_g( gO, BHD->g_ini[1], dgO2);

  VecSum(gH, &sumH);
  VecSum(gO, &sumO);

  VecGetArray(f, &f_vec);
  f_vec[0] = SQR(con*(sumH-sumO)-1.0);
  VecRestoreArray(f, &f_vec);

  return 0;
}


static real GetOptimalStepSize(State *BHD, Vec dg_newO, Vec dgO, Vec dgH)
{
  SNES snes;
  Vec s, f;
  PetscScalar *s_vec;
  StepData *SD;
  real step;

  SD = (StepData*) malloc(sizeof(*SD));
  SD->BHD = BHD;
  SD->dg_newO = dg_newO;
  SD->dgO = dgO;
  SD->dgH = dgH;

  VecCreateSeq(PETSC_COMM_SELF, 1, &s);
  VecCreateSeq(PETSC_COMM_SELF, 1, &f);

  SNESCreate(PETSC_COMM_WORLD, &snes);

  SNESSetFunction(snes, f, ComputeStepFunction, (void*)SD);
  SNESSetType(snes, SNESLS);

  SNESSolve(snes, PETSC_NULL, s);

  VecGetArray(s, &s_vec);
  step = s_vec[0];
  VecRestoreArray(s, &s_vec);

  VecDestroy(s);
  VecDestroy(f);
  SNESDestroy(snes);

  return step;
}

static void ComputedgFromg (Vec dg, Vec g0, Vec g)
{
  int local_size, i;
  PetscScalar *dg_vec, *g_vec, *g0_vec;
  real k;

  VecGetLocalSize(dg, &local_size);

  VecGetArray(dg, &dg_vec);
  VecGetArray(g0, &g0_vec);
  VecGetArray(g, &g_vec);
  for(i=0; i<local_size; i++)
    {
      k = g_vec[i];
      if(k<1.0e-8)
        k=1.0e-8;
      dg_vec[i] = -g0_vec[i] - log(k);

    }

  VecRestoreArray(dg, &dg_vec);
  VecRestoreArray(g0, &g0_vec);
  VecRestoreArray(g, &g_vec);
}

/*
 * Does the mixing:
 *
 *     dg := a * dg_new + (1 - a) * dg
 *
 * Returns the norm of the difference |dg_new - dg|.
 */
static real mix (Vec dg, Vec dg_new, real a, Vec work)
{
    real norm;

    /* Move dg */
    VecCopy(dg, work);

    /* dg' = a * dg_new + (1 - a) * dg */
    VecAXPBY(dg, a, (1-a), dg_new);

    /* work = dg - dg' = a * (dg - dg_new)

       That is why the divison by "a" below */
    VecAXPY(work, -1.0, dg);

    /* Norm of the change: */
    VecNorm(work, NORM_INFINITY, &norm);

    /* FIXME: why not computing |dg_new - dg| directly? Would be valid
       also for a == 0: */
    return norm / a;
}



/*
  This function is the main entry  point for the BGY3dM equation for a
  2-site solvent and an arbitrary solute.  The two vectors in

  Vec g[2], intent(out)

  are  initialized as global  distributed arrays  and filled  with the
  solvent site  distributions. It is the responsibility  of the caller
  to destroy them when no more needed.
 */
void bgy3d_solve_with_solute (const ProblemData *PD,
                              int n, const Site solute[n],
                              Vec g[2])
{
  real norm;
  Vec t_vec;                 /* used for all sites */
  Vec uc;                    /* Coulomb long, common for all sites. */
  Vec dg[2], dg_acc, work;
  PetscScalar dg_norm[2];
  int namecount = 0;
  char nameH[20], nameO[20];

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (2-site) equation ...\n");

  State BHD = initialize_state (PD);

  if (r_HH > 0.0)
      PetscPrintf(PETSC_COMM_WORLD,"WARNING: Solvent not a 2-Site model!\n");

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Mixing parameter: */
  const real a0 = PD->lambda;

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  /* Zeropad */
  const real zpad = PD->zpad;

  /* Inverse temperature: */
  const real beta = PD->beta;

  ImposeBoundaryCondition_Initialize (&BHD, zpad);
#ifdef L_BOUNDARY
  BHD.zpad = zpad;
  /* Assemble Laplacian matrix */
  InitializeLaplaceMatrix(&BHD, zpad);
  /* Create KSP environment */
  InitializeKSPSolver(&BHD);
#endif
#ifdef L_BOUNDARY_MG
  BHD.zpad = zpad;

  InitializeDMMGSolver(&BHD);
#endif

  /* These will hold  FFT of the current g. Allocate  enough to hold a
     local portion of the grid and free after the loop. */
  fftw_complex *g_fft[2];

  for (int i = 0; i < 2; i++) {
      g_fft[i] = bgy3d_fft_malloc (BHD.da);

      DACreateGlobalVector (BHD.da, &dg[i]);

      /* Here the storage for the output is allocated, the caller will
         have to destroy them: */
      DACreateGlobalVector (BHD.da, &g[i]);
  }

  /* These  are the  (four) kernels  HH, HO,  OH, OO.  Note that  HO =
     OH. */
  fftw_complex *ker_fft_S[2][2];
  fftw_complex *ker_fft_L[2][2];

  for (int i = 0; i < 2; i++)
      for (int j = 0; j <= i; j++) {
          ker_fft_S[i][j] = bgy3d_fft_malloc (BHD.da);
          ker_fft_S[j][i] = ker_fft_S[i][j];

          ker_fft_L[i][j] = bgy3d_fft_malloc (BHD.da);
          ker_fft_L[j][i] = ker_fft_L[i][j];
      }

  /* There is no point to  transform each contribution computed in the
     momentum scpace back, accumulate them on the k-grid here: */
  fftw_complex *dg_acc_fft = bgy3d_fft_malloc (BHD.da);

  DACreateGlobalVector(BHD.da, &dg_acc);
  DACreateGlobalVector(BHD.da, &work);
  DACreateGlobalVector(BHD.da, &t_vec); /* used for all sites */
  DACreateGlobalVector(BHD.da, &uc);    /* common for all sites */

  /* XXX: Here g0 = beta  * (VM_LJ + VM_coulomb_short) actually.  See:
          (5.106) and (5.108) in Jager's thesis. It is not filled with
          data yet, I assume. */
  Vec *g0 = BHD.g_ini;          /* FIXME: aliasing! */

  /* set initial guess*/
  VecSet(dg[0],0);
  VecSet(dg[1],0);

  if (bgy3d_getopt_test ("--from-g2")) {
      ComputedgFromg (dg[0], g0[0], BHD.g2[0][1]);
      ComputedgFromg (dg[1], g0[1], BHD.g2[1][1]);
  }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      dg[0] = bgy3d_load_vec ("dg0.bin"); /* dgH */
      dg[1] = bgy3d_load_vec ("dg1.bin"); /* dgO */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  }

  for (real damp = damp_start; damp <= 1; damp += 0.1) {

      /*
        FIXME: I guess the logic with damping factors can be made more
        straightforward:
      */
      real damp_LJ = (damp >= 0 ? 1.0 : 0.0); /* yes, >=, so the
                                                 original */

      /*
        Return F  * g2.  Note the  calculation of F is  divided due to
        long range Coulomb interation.   See comments in the function.
        Here F is force within solvents particles.

        In RecomputeInitialFFTs(), pairwise long-range interactions in
        BHD.u2[][] are computed by ComputeFFTfromCoulomb().

        According to Page. 115 - 116:

                  0
          g(x) = g (x) exp[-u(x)]

        with

           0               LJ       Cs       Cl
          g (x) = exp[-β (V  (x) + V  (x) + V  (x)]

        then g(x) rewritten as:

                 ~ 0            Cl
          g(x) = g  (x) exp[-β V  (x) - u(x)]

        with

          ~ 0               LJ       Cs
          g  (x) = exp[-β (V  (x) + V  (x))]

        then BGY equation written as:

           ~              Cl
          Δu = K(g) + β ΔV

        with
                 ~ 0            Cl              ~ 0         ~
          g(x) = g  (x) exp[-β V  (x) - u(x)] = g  (x) exp[-u(x)]

                         ~
        the solution  of u can  be repsented by a  difference of
        two functions:

          ~   -    *
          u = u - u

        while:

           -              Cl       -
          Δu = K(g) + β ΔV   in Ω, u(∂Ω) = f,

            *           *
          Δu = 0 in Ω, u (∂Ω) = f,

        so after solving:

          Δu = K(g)

        we get:

          -          Cl
          u = u + β V

         Cl         Cl
        V      and V      are needed to calculated beforehand
         (A, M)     (B, M)

        and  sum to  the solution  by the  end of  each  iteration for
        solvent  site A  and B,  in  the code  hereafter.  These  site
        specific potentials can be  obtained by multiplying the common
        long-range Coulomb potential (stored in Vec uc) by the solvent
        site charge.
      */

      RecomputeInitialFFTs(&BHD, (damp > 0.0 ? damp : 0.0), 1.0);

      /* FIXME:  avoid  storing  vectors  fg2XY* on  the  grid  across
         iterations. Only a scalar kernel is needed: */
      kernel (BHD.da, BHD.PD, BHD.f_g2_fft[0][0], NULL, ker_fft_S[0][0]);
      kernel (BHD.da, BHD.PD, BHD.f_g2_fft[0][1], NULL, ker_fft_S[0][1]); /* == [1][0] */
      kernel (BHD.da, BHD.PD, BHD.f_g2_fft[1][1], NULL, ker_fft_S[1][1]);

      kernel (BHD.da, BHD.PD, BHD.fl_g2_fft[0][0], BHD.u2_fft[0][0],  ker_fft_L[0][0]);
      kernel (BHD.da, BHD.PD, BHD.fl_g2_fft[0][1], BHD.u2_fft[0][1], ker_fft_L[0][1]); /* == [1][0] */
      kernel (BHD.da, BHD.PD, BHD.fl_g2_fft[1][1], BHD.u2_fft[1][1],  ker_fft_L[1][1]);

      /* FIXME: what is  the point to split the  kernel in two pieces?
         Redefine S := S + L and forget about L */
      for (int i = 0; i < 2; i++)
          for (int j = 0; j <= i; j++) {
              /* S := (damp / damp0) * L + damp_LJ * S */
              bgy3d_fft_axpby (BHD.da, ker_fft_S[i][j], damp / damp0, damp_LJ, ker_fft_L[i][j]);
          }

      /* Fill g0[0], g0[1] (alias BHD.g_ini[], also see the definition
         above) and  uc with VM_Coulomb_long.  No other  fields of the
         struct State except those passed explicitly are modified: */
      bgy3d_solute_field (&BHD,
                          2, solvent,
                          g0, uc, /* intent(out) */
                          n, solute,
                          (damp > 0.0 ? damp : 0.0), 1.0);

      /* Historically short-range  potential is scaled  by the inverse
         temperature: */
      for (int i = 0; i < 2; i++)
        VecScale (g0[i], beta);

      /* FIXME:  Check if this  is redundant  --- it  was mechanically
         moved from  the body of the  above func (because  it does not
         belong into solute code): */
      FOR_DIM {
        VecSet(BHD.F_l[0][0][dim], 0.0);
        VecSet(BHD.F_l[1][1][dim], 0.0);
        VecSet(BHD.F_l[0][1][dim], 0.0); /* What is it used for? */
      }

      PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);

      //Smooth_Function(&BHD, g0[1], zpad-1, zpad, 0.0);
      //Smooth_Function(&BHD, g0[0], zpad-1, zpad, 0.0);

      /* XXX: See  p116-177 in thesis:  Boundary Conditions  (5.107) -
             (5.110):  first  impose  boundary condistion  then  solve
             laplacian equation  and substrate from g0.   State BHD is
             not modified  by these calls. Note that  t_vec appears to
             be intent(out) in these calls, the value is ignored. */
      ImposeLaplaceBoundary (&BHD, g0[0], t_vec, BHD.x_lapl[0], zpad, NULL);
      ImposeLaplaceBoundary (&BHD, g0[1], t_vec, BHD.x_lapl[1], zpad, NULL);

      /* XXX: then correct g0 with boundary condition again. State BHD
              is not modified by these calls: */
      Zeropad_Function (&BHD, g0[1], zpad, 0.0);
      Zeropad_Function (&BHD, g0[0], zpad, 0.0);
      /* g=g0*exp(-dg) */

      ComputeH2O_g (g[0], g0[0], dg[0]);
      ComputeH2O_g (g[1], g0[1], dg[1]);

      real dgH_old = 0.0, dgO_old = 0.0; /* Not sure  if 0.0 as inital
                                            value is right.  */
      real a1 = a0;             /* loop-local variable */
      for (int iter = 0, mycount = 0, upwards = 0; iter < max_iter; iter++) {

          real a;               /* Used only inside the loop. */
          if (!(iter % 10) && iter > 0)
            a = a1;             /* This is taken  in iteration 10, 20,
                                   etc.  "a1" is  modified  during the
                                   loop. */
          else
            a = a0;             /* This  is taken  in  iterations 0-9,
                                   11-19, etc.  "a0" remains unchanged
                                   during the loop. */

          PetscPrintf(PETSC_COMM_WORLD,"iter %d: function norms: %e ", iter + 1, NORM_REG);

          /* The  functions  Compute_H2O_interS/_C() use  preallocated
             fftw_complex  arrays in  State BHD  for work  but  do not
             re-define  any  of  the  Vec(tors)  except  those  passed
             explicitly.

             Same    holds    for    Solve_NormalizationH2O_smallII(),
             ImposeLaplaceBoundary()  and  Zeropad_Function()  to  the
             best of my (limited) knowledge.

             The    factor   (damp/    damp0)   in    the    call   to
             Compute_H2O_interS_C()  looks  odd,  but  note  that  the
             Coulomb  field was  defined  having the  factor damp0  in
             RecomputeInitialFFTs(). FIXME: avoid this ugliness. */

          /* Compute FFT of g[] for all sites: */
          for (int i = 0; i < 2; i++)
              g_fft[i] =
                  ComputeFFTfromVec_fftw (BHD.da, BHD.fft_plan_fw,
                                          g[i], g_fft[i],
                                          BHD.fft_scratch); /* work array */

          /* for H, O in that order ... */
          for (int i = 0; i < 2; i++) {

              /* ... sum over H, O  in that order. LJ, short- and long
                 range Coulomb,  and a  so called strange  addition is
                 accounted   for   in    the   kernel.   First   clear
                 accumulator: */
              bgy3d_fft_set (BHD.da, dg_acc_fft, 0.0);

              for (int j = 0; j < 2; j++) {
                  /* This increments the accumulator: */
                  apply (BHD.da, ker_fft_S[i][j], g_fft[j], beta * BHD.rhos[j], dg_acc_fft);
              }

              /* Compute   IFFT   of   dg_acc_fft  for   the   current
                 site. Other contributions are added to the real space
                 dg_acc below: */
              ComputeVecfromFFT_fftw (BHD.da, BHD.fft_plan_bw,
                                      dg_acc, /* result */
                                      dg_acc_fft,
                                      BHD.fft_scratch);

              /* FIXME:   ugly  branch.    Very  specific   to  2-site
                 models. Literal  constants 0  and 1, how  comes?

                 Vec  t_vec,  is intent(out)  here.   Pass the  g_fft.
                 Earlier   version,  Solve_NormalizationH2O_smallII(),
                 did FFT itself wasting one FFT per site: */
              if (i == 0)
                  bgy3d_solve_normalization (&BHD, g_fft[i], r_HO, g[1], t_vec);
              else
                  bgy3d_solve_normalization (&BHD, g_fft[i], r_HO, g[0], t_vec);

              /* Vec t_vec is intent(in) and work is intent(out): */
              Compute_dg_H2O_intra_ln (&BHD, t_vec, r_HO, work);

              VecAXPY(dg_acc, 1.0, work);

              /* Add Coulomb  field uc scaled  by the site  chanrge to
                 the accumulator: */
              VecAXPY(dg_acc, solvent[i].charge, uc);

              /* Vec t_vec is intent(out) here: */
              ImposeLaplaceBoundary(&BHD, dg_acc, t_vec, BHD.x_lapl[i], zpad, NULL);
              Zeropad_Function(&BHD, dg_acc, zpad, 0.0);

              /*
               * Mix dg and dg_new with a fixed ration "a":
               *
               *     dg' = a * dg_new + (1 - a) * dg
               *     norm = |dg_new - dg|
               */
              dg_norm[i] = mix (dg[i], dg_acc, a, work); /* last arg is a temp */
          } /* over sites i */

          PetscPrintf(PETSC_COMM_WORLD,"H= %e (a=%f) ", dg_norm[0], a);
          PetscPrintf(PETSC_COMM_WORLD,"O= %e (a=%f) ", dg_norm[1], a);

          /* Now that dg[] has bee computed one can safely update g[].
             Compute g := exp[-(g0 + dg)], with a sanity check: */
          for (int i = 0; i < 2; i++)
              ComputeH2O_g (g[i], g0[i], dg[i]);

          /* Again a  strange case of  literal constants 0 and  1. Why
             not 1, 0? */
          norm = ComputeCharge(&BHD, g[0], g[1]);
          PetscPrintf(PETSC_COMM_WORLD, " %e ", norm);

          /*
           * Fancy step size control.
           *
           * FIXME: weired  logic.  The  code below appears  to modify
           * "upwards",  "a1", and  "mycount" by  eventually resetting
           * the latter to  zero. The goal might be  to tweak the real
           * coefficient  "a1"   depending  on  iteration   count  and
           * convergence.    Everytime  "mycount"   becomes   >20  the
           * coefficient "a1" is  changed. The compiler is complaining
           * that "upwards, dgH_old, dgO_old maybe used uninitialized"
           * here. At the moment I  am not able to confirm/reject that
           * claim.
           */
          mycount++;

          if (((iter - 1) % 10) &&
             (dgH_old < dg_norm[0] || dgO_old < dg_norm[1]))
            {
              upwards = 1;
            }
          else if (iter > 20 && !((iter - 1) % 10) && upwards == 0 &&
                  (dgH_old < dg_norm[0] || dgO_old < dg_norm[1]))
            {
              a1 /= 2.;
              if(a1 < a0)
                a1 = a0;
              mycount = 0;
            }
          else
            upwards = 0;

          if (mycount > 20) {
            /*
             * Scale  the coefficient "a1"  up by  a factor,  but make
             * sure it is not above 1.0. Reset mycount.
             */
            if (a1 <= 0.5) {
              a1 *= 2.0;
            }
            else {
              a1 = 1.0;
            }
            mycount = 0;
          }
          /* otherwise leave "a1" and "mycount" unchanged */

          PetscPrintf(PETSC_COMM_WORLD,"count= %d  upwards= %d", mycount, upwards);
          dgH_old = dg_norm[0];
          dgO_old = dg_norm[1];

          /*********************************/

          PetscPrintf(PETSC_COMM_WORLD,"\n");

          if (dg_norm[0] <= norm_tol &&  dg_norm[1] <= norm_tol) //&& NORM_REG<5.0e-2)
            break;

      } /* iter loop */
      /*************************************/

      /* FIXME:  Debug  output  from  every iteration  with  different
         overall  scale  factors  damp/damp_LJ.  Remove when  no  more
         needed. */
      namecount++;
      sprintf(nameH, "vec0-%d.m", namecount-1);
      sprintf(nameO, "vec1-%d.m", namecount-1);

      PetscPrintf(PETSC_COMM_WORLD,"Writing files...");
      bgy3d_save_vec_ascii (nameH, g[0]); /* g_H */
      bgy3d_save_vec_ascii (nameO, g[1]); /* g_O */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
      /************************************/

      /* Save dg to binary file. FIXME: Why dg and not g? */
      if (bgy3d_getopt_test ("--save-H2O")) {
          PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
          bgy3d_save_vec ("dg0.bin", dg[0]); /* dgH */
          bgy3d_save_vec ("dg1.bin", dg[1]); /* dgO */
          PetscPrintf(PETSC_COMM_WORLD,"done.\n");
      }
  } /* damp loop */

  /* Clean up and exit ... */
  bgy3d_fft_free (dg_acc_fft);

  for (int i = 0; i < 2; i++) {
      /* Delegated to the caller:
         VecDestroy (g[i]); */
      VecDestroy (dg[i]);

      bgy3d_fft_free (g_fft[i]);

      for (int j = 0; j <= i; j++) {
          bgy3d_fft_free (ker_fft_S[i][j]);
          bgy3d_fft_free (ker_fft_L[i][j]);
      }
  }

  VecDestroy (dg_acc);
  VecDestroy (work);
  VecDestroy (t_vec);
  VecDestroy (uc);

  finalize_state (&BHD);
}

/* This one emulates historical solver interface: */
Vec BGY3dM_solve_H2O_2site(const ProblemData *PD, Vec g_ini)
{
  (void) g_ini;                 /* FIXME: interface obligation */

  int n;                        /* number of solute sites */
  const Site *sites;            /* [n], array of sites */
  const char *name;             /* human readable name */
  Vec g[2];                     /* solution */

  /* Solute index between 0 and 5 (inclusive): */
  const int solute = PD->solute;

  /* Get the solute from the tables: */
  bgy3d_solute_get (solute, &n, &sites, &name);

  /* Code used to be verbose: */
  PetscPrintf(PETSC_COMM_WORLD,"Solute is %s.\n", name);

  /* This does the  real work. Vec g[2] is  intent(out) in all senses,
     dont forget to destroy them: */
  bgy3d_solve_with_solute (PD, n, sites, g);

  /* Save final distribution, use binary format: */
  bgy3d_save_vec ("g0.bin", g[0]); /* gH */
  bgy3d_save_vec ("g1.bin", g[1]); /* gO */

  for (int i = 0; i < 2; i++)
    VecDestroy (g[i]);

  return PETSC_NULL;            /* fake, interface obligation */
}

Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini)
{
  (void) g_ini;                 /* FIXME: interface obligation */

  real a1, a, damp, damp_LJ;
  real count=0.0, norm, aO;
  int iter;
  Vec g0H, g0O, dgH, dgO,  dg_new, dg_new2, f, gH, gO;
  Vec tH, tO, dg_newH, dg_newO;
  Vec uc;                       /* common for all sites */
  PetscScalar dgH_norm, dgO_norm;
  real dgH_old, dgO_old;
  int mycount=0, upwards, namecount=0;
  char nameH[20], nameO[20];
  real ti;
  int iteri;

  Vec dg_histO, dg_histH;

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (3-site) equation ...\n");

  State BHD = initialize_state (PD);
  BHD.rhos[0] = 2.*BHD.rhos[0];

  if (r_HH < 0.0) {
    PetscPrintf(PETSC_COMM_WORLD,"Solvent not a 3-Site model!\n");
    exit(1);
  }

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Mixing parameter: */
  real a0 = PD->lambda;         /* not const */

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  /* Zeropad */
  const real zpad = PD->zpad;

  ImposeBoundaryCondition_Initialize( &BHD, zpad);
#ifdef L_BOUNDARY
  BHD.zpad = zpad;
  /* Assemble Laplacian matrix */
  InitializeLaplaceMatrix(&BHD, zpad);
  /* Create KSP environment */
  InitializeKSPSolver(&BHD);
#endif
#ifdef L_BOUNDARY_MG
  BHD.zpad = zpad;
  /* Create KSP environment */
  InitializeDMMGSolver(&BHD);
#endif

  DACreateGlobalVector(BHD.da, &gH);
  DACreateGlobalVector(BHD.da, &gO);
  DACreateGlobalVector(BHD.da, &dgH);
  DACreateGlobalVector(BHD.da, &dgO);
  DACreateGlobalVector(BHD.da, &dg_new);
  DACreateGlobalVector(BHD.da, &dg_new2);
  DACreateGlobalVector(BHD.da, &f);

  DACreateGlobalVector(BHD.da, &tH);
  DACreateGlobalVector(BHD.da, &tO);

  DACreateGlobalVector(BHD.da, &dg_newH);
  DACreateGlobalVector(BHD.da, &dg_newO);

  DACreateGlobalVector(BHD.da, &uc); /* common for all sites */

  DACreateGlobalVector(BHD.da, &dg_histO);
  DACreateGlobalVector(BHD.da, &dg_histH);
  VecSet(dg_histH, 0.0);
  VecSet(dg_histO, 0.0);

  g0H=BHD.g_ini[0];
  g0O=BHD.g_ini[1];

/*   VecView(BHD.F_l[0][1][0],PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dg_new,0.0);

  if (bgy3d_getopt_test ("--from-g2")) {
      ComputedgFromg (dgH, g0H, BHD.g2[0][1]);
      ComputedgFromg (dgO, g0O, BHD.g2[1][1]);
  }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      dgH = bgy3d_load_vec ("dg0.bin"); /* dgH */
      dgO = bgy3d_load_vec ("dg1.bin"); /* dgO */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  }

  for( damp=damp_start; damp <=1; damp+=0.1)
    {
      if (damp <= -0.01)
        {
          damp_LJ = 0.0;
        }
      else if (damp==0.0)
        {
          damp_LJ = 1.0;
        }
      else
        {
          damp_LJ = 1.0;
          count += 1.0;
          a0 = 0.1 / (count + 5.0);
        }

      RecomputeInitialFFTs(&BHD, (damp > 0.0 ? damp : 0.0), damp_LJ);

      /* Compute solute field in this block:*/
      {
        int n;                        /* number of solute sites */
        const Site *sites;            /* [n], array of sites */
        const char *name;             /* human readable name */

        /* Get the solute from the tables: */
        bgy3d_solute_get (4, &n, &sites, &name); /* Butanoic Acid */

        /* This does the real work: */
        bgy3d_solute_field (&BHD,
                            2, solvent,
                            BHD.g_ini, uc, /* intent(out) */
                            n, sites,
                            (damp > 0.0 ? damp : 0.0), damp_LJ);
      }

      PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);

      /* Historically   short-range  potential   is   stored  with   a
         factor: */
      for (int i = 0; i < 2; i++)
        VecScale (BHD.g_ini[i], BHD.PD->beta);

      /* FIXME:  Check if this  is redundant  --- it  was mechanically
         moved from the body of the bgy3d_solute_field() func (because
         it does not belong into solute code): */
      FOR_DIM {
        VecSet(BHD.F_l[0][0][dim], 0.0);
        VecSet(BHD.F_l[1][1][dim], 0.0);
        VecSet(BHD.F_l[0][1][dim], 0.0); /* What is it used for? */
      }

      //Smooth_Function(&BHD, g0O, zpad-1, zpad, 0.0);
      //Smooth_Function(&BHD, g0H, zpad-1, zpad, 0.0);

      ImposeLaplaceBoundary(&BHD, g0H, tH, BHD.x_lapl[0], zpad, NULL);
      ImposeLaplaceBoundary(&BHD, g0O, tH, BHD.x_lapl[1], zpad, NULL);
      Zeropad_Function(&BHD, g0O, zpad, 0.0);
      Zeropad_Function(&BHD, g0H, zpad, 0.0);
      /* g=g0*exp(-dg) */

      ComputeH2O_g( gH, g0H, dgH);
      ComputeH2O_g( gO, g0O, dgO);


  /*     VecCopy(BHD.g2[0][1], gH);  */
/*       VecCopy(BHD.g2[1][1], gO);  */

      a=a0;
      a1=a0;
      for(iter=0; iter<max_iter; iter++)
        {

          if( !(iter%10) && iter>0 )
            a=a1;
          else
            a=a0;

/*        if( !(iter%50) && iter>0 && NORM_REG>=5.0e-2) */
/*          NORM_REG/=2.; */

          PetscPrintf(PETSC_COMM_WORLD,"iter %d: function norms: %e ", iter+1, NORM_REG);


          /* H */
          VecSet(dg_new,0.0);
          Compute_H2O_interS(&BHD, BHD.f_g2_fft[0][1], gO, BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.f_g2_fft[0][0], gH, BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][1], gO, BHD.u2_fft[0][1], (damp / damp0) * BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][0], gH, BHD.u2_fft[0][0], (damp / damp0) * BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_smallII( &BHD, gH, r_HH, gH, tH , dg_new2, f, zpad);
          //Compute_dg_H2O_intra_lnIII(&BHD, gH, tH, r_HH, dg_new2, f);
          //Compute_dg_H2O_intra_lnII(&BHD, gH, tH, r_HH, dg_new2, f);
          Compute_dg_H2O_intra_ln(&BHD, tH, r_HH, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);
          Solve_NormalizationH2O_smallII( &BHD, gH, r_HO, gO, tO , dg_new2, f, zpad);
          //Compute_dg_H2O_intra_lnIII(&BHD, gO, tO, r_HO, dg_new2, f);
          //Compute_dg_H2O_intra_lnII(&BHD, gO, tO, r_HO, dg_new2, f);
          Compute_dg_H2O_intra_ln(&BHD, tO, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

/*        VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);      */
/*        exit(1);    */

          //dgH_norm = ComputeCharge(&BHD, gH, gO);
          //PetscPrintf(PETSC_COMM_WORLD, " %e ", dgH_norm);
/*        if(dgH_norm<1.0) */
/*          VecScale(dg_new, 1.2); */
/*        else if(dgH_norm >1.0) */
/*          VecScale(dg_new, 0.8); */

          /* Copy  the electrostatic potential  scaled by  the solvent
             site charges into predefined locations: */
          VecAXPY(dg_new, solvent[0].charge, uc);
/*        dgH_norm = ImposeBoundaryConditionII(&BHD, dg_new, zpad); */
/*        PetscPrintf(PETSC_COMM_WORLD, " %e ", dgH_norm); */
          //VecShift(dg_new, -dgH_norm);
          //ImposeBoundaryCondition( &BHD, dg_new);

          ti=ImposeLaplaceBoundary(&BHD, dg_new, tH, BHD.x_lapl[0], zpad, &iteri);
          Zeropad_Function(&BHD, dg_new, zpad, 0.0);
          //Smooth_Function(&BHD, dg_new, zpad-1, zpad, 0.0);

          PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);



          VecCopy(dg_new, dg_newH);

/*        VecCopy(dgH, f); */
/*        VecAXPBY(dgH, a, (1-a), dg_new); */
/*        VecAXPY(f, -1.0, dgH); */
/*        VecNorm(f, NORM_INFINITY, &dgH_norm); */
/*        PetscPrintf(PETSC_COMM_WORLD,"H= %e  ", dgH_norm/a); */
          //ComputeH2O_g( gH, g0H, dgH);


          /* O */
          VecSet(dg_new,0.0);
          Compute_H2O_interS(&BHD, BHD.f_g2_fft[1][1], gO, BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.f_g2_fft[0][1], gH, BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[1][1], gO, BHD.u2_fft[1][1], (damp / damp0) * BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][1], gH, BHD.u2_fft[0][1], (damp / damp0) * BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);



          Solve_NormalizationH2O_smallII( &BHD, gO, r_HO, gH, tH , dg_new2, f, zpad);
          Compute_dg_H2O_intra_ln(&BHD, tH, r_HO, dg_new2);
          VecCopy (dg_new2, f); /* FIXME: need that? */
          //Compute_dg_H2O_intra_lnII(&BHD, gH, tH, r_HO, dg_new2, f);
          //Compute_dg_H2O_intra_lnIII(&BHD, gH, tH, r_HO, dg_new2, f);
          VecAXPY(dg_new, 2.0, dg_new2);

/*        VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);      */
/*        exit(1);    */


          //dgO_norm = ComputeCharge(&BHD, gH, gO);
          //PetscPrintf(PETSC_COMM_WORLD, " %e ", dgO_norm);
/*        if(dgO_norm<1.0) */
/*          VecScale(dg_new, 0.9); */
/*        else if(dgO_norm >1.0) */
/*          VecScale(dg_new, 1.1); */

          /* Copy  the electrostatic potential  scaled by  the solvent
             site charges into predefined locations: */
          VecAXPY(dg_new, solvent[1].charge, uc);
/*        dgO_norm = ImposeBoundaryConditionII(&BHD, dg_new, zpad); */
/*        PetscPrintf(PETSC_COMM_WORLD, " %e ", dgO_norm); */
          //VecShift(dg_new, -dgO_norm);
          //ImposeBoundaryCondition( &BHD, dg_new);

          //ti=ImposeLaplaceBoundary(&BHD, dg_new, tH, tO, zpad, &iteri);

          ti=ImposeLaplaceBoundary(&BHD, dg_new, tH, BHD.x_lapl[1], zpad, &iteri);
          Zeropad_Function(&BHD, dg_new, zpad, 0.0);
          //Smooth_Function(&BHD, dg_new, zpad-1, zpad, 0.0);
          PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);

/*        if(iter==51) */
/*          { */
/*            VecView(tO,PETSC_VIEWER_STDERR_WORLD);  */
/*            exit(1); */
/*          }  */

          VecCopy(dg_new, dg_newO);


/*        VecCopy(dgO, f); */
/*        VecAXPBY(dgO, a, (1-a), dg_new); */
/*        VecAXPY(f, -1.0,  dgO); */
/*        VecNorm(f, NORM_INFINITY, &dgO_norm); */
/*        PetscPrintf(PETSC_COMM_WORLD,"O= %e  ", dgO_norm/a); */
          //ComputeH2O_g( gO, g0O, dgO);

/*        if(iter<=10) */
/*          a=0.01; */
/*        else */
/*          a=0.05; */
/*        if( iter==0) */
/*          { */
/*            VecScale(dg_newO, a); */
/*            VecScale(dg_newH, a); */
/*            EnforceNormalizationCondition(&BHD, dg_newO, dg_newH, gO, gH); */
/*            a=1.0; */
/*          } */

          /* Move dgH */
          VecCopy(dgH, f);
/*        VecAXPY(dg_histH, 1.0, dg_newH); */
/*        norm = (iter)%50 + 1; */

          VecAXPBY(dgH, a, (1-a), dg_newH);
          VecAXPY(f, -1.0, dgH);
          VecNorm(f, NORM_INFINITY, &dgH_norm);

          PetscPrintf(PETSC_COMM_WORLD,"H= %e (a=%f) ", dgH_norm/a, a);

          /* Move dgO */
          if(0&& iter>10)
            {
              aO = GetOptimalStepSize(&BHD, dg_newO, dgO, dgH);
              //if(aO<0) aO = a;
              VecCopy(dgO, f);
              VecAXPBY(dgO, SQR(aO), (1-SQR(aO)), dg_newO);
              //VecCopy(dg_newO, dgO);
              //VecScale(dgO, aO);

              VecAXPY(f, -1.0,  dgO);
              VecNorm(f, NORM_INFINITY, &dgO_norm);
              PetscPrintf(PETSC_COMM_WORLD,"O= %e (a=%f) ", dgO_norm/aO, aO);
            }
          else
            {
              VecCopy(dgO, f);
/*            VecAXPY(dg_histO, 1.0, dg_newO); */
/*            norm = (iter)%20 + 1; */
              VecAXPBY(dgO, a, (1-a), dg_newO);
              VecAXPY(f, -1.0,  dgO);
              VecNorm(f, NORM_INFINITY, &dgO_norm);
              PetscPrintf(PETSC_COMM_WORLD,"O= %e (a=%f) ", dgO_norm/a, a);
            }
          //a=a+0.002;
/*        if( !((iter+1)%50) && iter>0) */
/*          { */
/*            a=0.01; */
/*            PetscPrintf(PETSC_COMM_WORLD,"Adding history_vec..."); */
/*            VecAXPY(dgO, -2.0, dg_histO); */
/*            VecAXPY(dgH, -2.0, dg_histH); */
/*            VecSet(dg_histO, 0.0); */
/*            VecSet(dg_histH, 0.0); */
/*          } */
          ComputeH2O_g( gH, g0H, dgH);
          ComputeH2O_g( gO, g0O, dgO);
          norm = ComputeCharge(&BHD, gH, gO);
          PetscPrintf(PETSC_COMM_WORLD, " %e ", norm);
          //EnforceNormalizationCondition(&BHD, dgO, dgH, gO, gH);


/*        if( !(iter%10) &&iter>0  ) */
/*          { */
/*            EnforceNormalizationCondition(&BHD, dgO, dgH, gO, gH); */

/*          } */

/*        if(dgO_norm/a<0.1 && dgH_norm/a<0.1) */
/*          { */
/*            EnforceNormalizationCondition(&BHD, dgO, dgH, gO, gH); */
/*            a0=0.1; */
/*          } */

/*        if( fabs(norm-1.0) > 0.05 ) */
/*          a=0.005; */
/*        else  */
/*          a=0.01; */
           /* (fancy) step size control */
          mycount++;
          if( ((iter-1)%10) &&
              (dgH_old<dgH_norm/a || dgO_old<dgO_norm/a ) )
            {
              upwards=1;
            }
          else if(iter>20 && !((iter-1)%10) && upwards==0 &&
                  (dgH_old<dgH_norm/a || dgO_old<dgO_norm/a ) )
            {
              a1 /=2.;
              if(a1<a0)
                a1=a0;
              mycount=0;
            }
          else
            upwards=0;

          if(mycount>20 && a1<=0.5)
            {
              a1*=2.;
              mycount=0;
            }
          else if(mycount>20 && a1>0.5)
            {
              a1=1.0;
              mycount=0;
            }
          PetscPrintf(PETSC_COMM_WORLD,"count= %d  upwards= %d", mycount, upwards);
          dgH_old = dgH_norm/a;
          dgO_old = dgO_norm/a;

          /*********************************/

          PetscPrintf(PETSC_COMM_WORLD,"\n");

          if(dgH_norm/a<=norm_tol &&  dgO_norm/a<=norm_tol) //&& NORM_REG<5.0e-2)
            break;

        }
      /*************************************/
      /* output */
      namecount++;
      sprintf(nameH, "vec0-%d.m", namecount-1);
      sprintf(nameO, "vec1-%d.m", namecount-1);

      PetscPrintf(PETSC_COMM_WORLD,"Writing files...");
      bgy3d_save_vec_ascii (nameH, gH); /* g_H */
      bgy3d_save_vec_ascii (nameO, gO); /* g_O */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");

      /************************************/

      /* save g to binary file */
      if (bgy3d_getopt_test ("--save-H2O")) {
          PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
          bgy3d_save_vec ("dg0.bin", dgH); /* dgH */
          bgy3d_save_vec ("dg1.bin", dgH); /* dgO */
          PetscPrintf(PETSC_COMM_WORLD,"done.\n");
      }
    }




  VecDestroy(gH);
  VecDestroy(gO);
  VecDestroy(dgH);
  VecDestroy(dgO);
  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(f);

  VecDestroy(tH);
  VecDestroy(tO);

  VecDestroy (uc);

  VecDestroy(dg_newH);
  VecDestroy(dg_newO);
  VecDestroy(dg_histH);
  VecDestroy(dg_histO);

  finalize_state (&BHD);

  return PETSC_NULL;
}
