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
#include "bgy3d_MG.h"           /* InitializeDMMGSolver */
#endif

#ifdef WITH_COMPLEX
#include <complex.h>
#endif

#define damp0 0.001

extern real NORM_REG;

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
  bgy3d_fft_init_da (PD->N, &(BHD.fft_plan_fw), &(BHD.fft_plan_bw), &(BHD.da), NULL);
#else
  /* multigrid, apparently needs two descriptors: */
  bgy3d_fft_init_da (PD->N, &(BHD.fft_plan_fw), &(BHD.fft_plan_bw), &(BHD.da), &(BHD.da_dmmg));
#endif

  const DA da = BHD.da;         /* shorter alias */

  /* Create global vectors */
  ierr = DACreateGlobalVector(da, &(BHD.g_ini[0])); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.g_ini[1])); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.gHO_ini)); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.uc[0])); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.uc[1])); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.ucHO)); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.g2H)); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.g2O)); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.g2HO)); assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD.pre)); assert (!ierr);

  FOR_DIM
    {
      ierr = DACreateGlobalVector(da, &(BHD.fH[dim])); assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD.fO[dim])); assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD.fHO[dim])); assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD.fH_l[dim])); assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD.fO_l[dim])); assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD.fHO_l[dim])); assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD.v[dim])); assert (!ierr);
    }


  VecSet(BHD.pre,0.0);

#ifdef L_BOUNDARY
   /* Create Matrix with appropriate non-zero structure */
  DAGetMatrix( da, MATMPIAIJ, &(BHD.M));
  DACreateGlobalVector(da, &(BHD.x_lapl[0]));
  DACreateGlobalVector(da, &(BHD.x_lapl[1]));
  VecSet(BHD.x_lapl[0], 0.0);
  VecSet(BHD.x_lapl[1], 0.0);
#endif

  /* Allocate memory for fft */
  FOR_DIM
    {
      BHD.fg2HH_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fg2OO_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fg2HO_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fg2HHl_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fg2OOl_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fg2HOl_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fg2_fft[dim] = bgy3d_fft_malloc (da);

      BHD.fO_fft[dim] = bgy3d_fft_malloc (da);
      BHD.fH_fft[dim] = bgy3d_fft_malloc (da);
    }

  BHD.g_fft = bgy3d_fft_malloc (da);
  BHD.gfg2_fft = bgy3d_fft_malloc (da);
  BHD.fft_scratch = bgy3d_fft_malloc (da);
  BHD.ucH_fft = bgy3d_fft_malloc (da);
  BHD.ucO_fft = bgy3d_fft_malloc (da);
  BHD.ucHO_fft = bgy3d_fft_malloc (da);
  BHD.wHO_fft = bgy3d_fft_malloc (da);
  BHD.wHH_fft = bgy3d_fft_malloc (da);



  /* Read g^2  from file */
/*   ReadPairDistribution(&BHD, "g2_OO", BHD.g2O); */
/*   ReadPairDistribution(&BHD, "g2_HH", BHD.g2H); */
/*   ReadPairDistribution(&BHD, "g2_HO", BHD.g2HO); */

#ifdef CS2
  ReadPairDistribution(&BHD, "g2C", BHD.g2O);
  ReadPairDistribution(&BHD, "g2S", BHD.g2H);
  ReadPairDistribution(&BHD, "g2CS", BHD.g2HO);
#else
  bgy3d_load_vec ("g00.bin", &BHD.g2H);
  bgy3d_load_vec ("g11.bin", &BHD.g2O);
  bgy3d_load_vec ("g01.bin", &BHD.g2HO);
#endif

  return BHD;
}



static void finalize_state (State *BHD)
{
  MPI_Barrier( PETSC_COMM_WORLD);

  FOR_DIM
    {
      VecDestroy(BHD->fH[dim]);
      VecDestroy(BHD->fO[dim]);
      VecDestroy(BHD->fHO[dim]);
      VecDestroy(BHD->fH_l[dim]);
      VecDestroy(BHD->fO_l[dim]);
      VecDestroy(BHD->fHO_l[dim]);
      VecDestroy(BHD->v[dim]);
      bgy3d_fft_free (BHD->fg2HH_fft[dim]);
      bgy3d_fft_free (BHD->fg2OO_fft[dim]);
      bgy3d_fft_free (BHD->fg2HO_fft[dim]);
      bgy3d_fft_free (BHD->fg2HHl_fft[dim]);
      bgy3d_fft_free (BHD->fg2OOl_fft[dim]);
      bgy3d_fft_free (BHD->fg2HOl_fft[dim]);
      bgy3d_fft_free (BHD->fg2_fft[dim]);

      bgy3d_fft_free (BHD->fO_fft[dim]);
      bgy3d_fft_free (BHD->fH_fft[dim]);
    }
  bgy3d_fft_free (BHD->g_fft);
  bgy3d_fft_free (BHD->gfg2_fft);
  bgy3d_fft_free (BHD->fft_scratch);
  bgy3d_fft_free (BHD->ucH_fft);
  bgy3d_fft_free (BHD->ucO_fft);
  bgy3d_fft_free (BHD->ucHO_fft);
  bgy3d_fft_free (BHD->wHO_fft);
  bgy3d_fft_free (BHD->wHH_fft);

  VecDestroy(BHD->g_ini[0]);
  VecDestroy(BHD->g_ini[1]);
  VecDestroy(BHD->gHO_ini);
  VecDestroy(BHD->uc[0]);
  VecDestroy(BHD->uc[1]);
  VecDestroy(BHD->ucHO);
  VecDestroy(BHD->g2H);
  VecDestroy(BHD->g2O);
  VecDestroy(BHD->g2HO);
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
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

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
  KSPCreate( PETSC_COMM_WORLD, &(BHD->ksp));
  KSPGetPC(BHD->ksp, &(pc));
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
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

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

  while( fscanf(fp,"%lf %lf", &(xg[index]), &(g[index])) == 2)
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
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

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

void RecomputeInitialFFTs(State *BHD, real damp, real damp_LJ)
{
  DA da;
  PetscScalar ***(fH_vec[3]),***(fO_vec[3]),***(fHO_vec[3]);
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]),***(fHOl_vec[3]);
  // PetscScalar ***gO_vec, ***gH_vec, ***gHO_vec;
  PetscScalar ***wHO_vec, ***wHH_vec;
  real r[3], r_s, h[3], interval[2], wconst_HO, wconst_HH, wG;
  int x[3], n[3], i[3];
  real epsilonH, epsilonO, epsilonHO;
  real sigmaH, sigmaO, sigmaHO;
  real q2H, q2O, q2HO;

  epsilonH = BHD->LJ_paramsH[0];
  epsilonO = BHD->LJ_paramsO[0];
  epsilonHO = BHD->LJ_paramsHO[0];

  sigmaH = BHD->LJ_paramsH[1];
  sigmaO = BHD->LJ_paramsO[1];
  sigmaHO = BHD->LJ_paramsHO[1];

  q2H = BHD->LJ_paramsH[2];
  q2O = BHD->LJ_paramsO[2];
  q2HO = BHD->LJ_paramsHO[2];


  const ProblemData *PD = BHD->PD;
  da = BHD->da;

  PetscPrintf(PETSC_COMM_WORLD,"Recomputing FFT data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);


  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];

  wG = 1./h[0];
  wconst_HO  =
    4.*M_PI*2.*(sqrt(M_PI)/4./pow(wG,3)-r_HO/2./SQR(wG)+SQR(r_HO)*sqrt(M_PI)/2./wG);
  wconst_HO = 1./wconst_HO;
  wconst_HH  =
    4.*M_PI*2.*(sqrt(M_PI)/4./pow(wG,3)-r_HH/2./SQR(wG)+SQR(r_HH)*sqrt(M_PI)/2./wG);
  wconst_HH = 1./wconst_HH;


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  FOR_DIM
    {
      VecSet(BHD->fH[dim],0.0);
      VecSet(BHD->fO[dim],0.0);
      VecSet(BHD->fHO[dim],0.0);
      VecSet(BHD->fH_l[dim],0.0);
      VecSet(BHD->fO_l[dim],0.0);
      VecSet(BHD->fHO_l[dim],0.0);
    }

  /* Compute Coulomb from fft part */
/*   ComputeFFTfromCoulombII(BHD, BHD->fO, BHD->fO_l, BHD->ucO_fft, BHD->LJ_paramsO, damp); */
/*   ComputeFFTfromCoulombII(BHD, BHD->fH, BHD->fH_l, BHD->ucH_fft, BHD->LJ_paramsH, damp); */
/*   ComputeFFTfromCoulombII(BHD, BHD->fHO, BHD->fHO_l, BHD->ucHO_fft, BHD->LJ_paramsHO, damp); */
/* XXX: uc[1] and uc[0] will be updated by RecomputeInitialSoluteData_XXX(),
 *      ucHO is not used at all */
  ComputeFFTfromCoulomb(BHD, BHD->uc[1], BHD->fO_l, BHD->ucO_fft, q2O, damp0);
  ComputeFFTfromCoulomb(BHD, BHD->uc[0], BHD->fH_l, BHD->ucH_fft, q2H, damp0);
  ComputeFFTfromCoulomb(BHD, BHD->ucHO, BHD->fHO_l, BHD->ucHO_fft, q2HO, damp0);
/*   FOR_DIM */
/*     { */
/*       VecAXPY(BHD->fO[dim], 1.0, BHD->fO_l[dim]); */
/*       VecAXPY(BHD->fH[dim], 1.0, BHD->fH_l[dim]); */
/*       VecAXPY(BHD->fHO[dim], 1.0, BHD->fHO_l[dim]); */
/*     } */

  FOR_DIM
    {
      DAVecGetArray(da, BHD->fH[dim], &(fH_vec[dim]));
      DAVecGetArray(da, BHD->fO[dim], &(fO_vec[dim]));
      DAVecGetArray(da, BHD->fHO[dim], &(fHO_vec[dim]));
      DAVecGetArray(da, BHD->fH_l[dim], &(fHl_vec[dim]));
      DAVecGetArray(da, BHD->fO_l[dim], &(fOl_vec[dim]));
      DAVecGetArray(da, BHD->fHO_l[dim], &(fHOl_vec[dim]));
    }
/*   VecSet(BHD->g2H,0.0); */
/*   VecSet(BHD->g2O,0.0); */
/*   VecSet(BHD->g2HO,0.0); */
/*   DAVecGetArray(da, BHD->g2H, &gH_vec); */
/*   DAVecGetArray(da, BHD->g2O, &gO_vec); */
/*   DAVecGetArray(da, BHD->g2HO, &gHO_vec); */

  DAVecGetArray(da, BHD->v[0], &wHO_vec);
  DAVecGetArray(da, BHD->v[1], &wHH_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          /* set force vectors */

          FOR_DIM
            r[dim] = i[dim]*h[dim]+interval[0];


          r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

            /* Lennard-Jones */
/*        gH_vec[i[2]][i[1]][i[0]] +=  */
/*          exp(-damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsH)); */
/*        gO_vec[i[2]][i[1]][i[0]] +=  */
/*          exp(-damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsO)); */
/*        gHO_vec[i[2]][i[1]][i[0]] +=  */
/*          exp(-damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsHO)); */


          /* omega_HH + omega_HO */
/*        wHO_vec[i[2]][i[1]][i[0]] = wconst_HO * exp(-SQR(wG)*SQR(r_s-r_HO)); */
/*        wHH_vec[i[2]][i[1]][i[0]] = wconst_HH * exp(-SQR(wG)*SQR(r_s-r_HH)); */
/*        wHO_vec[i[2]][i[1]][i[0]] = damp * Coulomb_long( r_s, BHD->LJ_paramsO); */

           FOR_DIM
            {
              /* Lennard-Jones */
              fH_vec[dim][i[2]][i[1]][i[0]] +=
                damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonH, sigmaH);
              fO_vec[dim][i[2]][i[1]][i[0]] +=
                damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonO, sigmaO);
              fHO_vec[dim][i[2]][i[1]][i[0]] +=
                damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonHO, sigmaHO);



              /* Coulomb short */
              fH_vec[dim][i[2]][i[1]][i[0]] +=
                damp * Coulomb_short_grad( r_s, r[dim], q2H);
              fO_vec[dim][i[2]][i[1]][i[0]] +=
                damp * Coulomb_short_grad( r_s, r[dim], q2O);
              fHO_vec[dim][i[2]][i[1]][i[0]] +=
                damp * Coulomb_short_grad( r_s, r[dim], q2HO);

              /* Coulomb long */
/*            fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/*            fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/*            fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

              /* Coulomb long */
 /*           fHl_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/*            fOl_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/*            fHOl_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

              /* Coulomb */
/*            fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsH); */
/*            fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsO); */
/*            fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/*              damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsHO); */


            }
        }

  FOR_DIM
    {
      DAVecRestoreArray(da, BHD->fH[dim], &(fH_vec[dim]));
      DAVecRestoreArray(da, BHD->fO[dim], &(fO_vec[dim]));
      DAVecRestoreArray(da, BHD->fHO[dim], &(fHO_vec[dim]));
      DAVecRestoreArray(da, BHD->fH_l[dim], &(fHl_vec[dim]));
      DAVecRestoreArray(da, BHD->fO_l[dim], &(fOl_vec[dim]));
      DAVecRestoreArray(da, BHD->fHO_l[dim], &(fHOl_vec[dim]));
    }
/*   DAVecRestoreArray(da, BHD->g2H, &gH_vec); */
/*   DAVecRestoreArray(da, BHD->g2O, &gO_vec); */
/*   DAVecRestoreArray(da, BHD->g2HO, &gHO_vec); */
  DAVecRestoreArray(da, BHD->v[0], &wHO_vec);
  DAVecRestoreArray(da, BHD->v[1], &wHH_vec);

/*   VecView(BHD->v[0],PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */

  /* Compute FFT(F*g^2) */
  // XXX: F*g2 = (F_LJ + F_coulomb_short) * g2 + (F_coulomb_long * g2 - F_coulomb_long) + F_coulomb_long
  // XXX: see (5.101) and (5.102) in Jager's thesis
  // XXX: FFT(F_coulomb_long) has been calculated as BHD->ucH_fft, BHD->ucO_fft, BHD->ucHO_fft by ComputeFFTfromCoulomb() above
  FOR_DIM
    {
      /* OO */
      // XXX: (F_LJ + F_coulomb_short) * g2 
      VecPointwiseMult(BHD->v[dim], BHD->g2O, BHD->fO[dim]);
      //VecAXPY(BHD->v[dim], -1.0, BHD->fO_l[dim]);
      // XXX: FFT((F_LJ + F_coulomb_short) * g2) 
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2OO_fft[dim], BHD->fft_scratch);
      /* Coulomb long */
      // XXX: F_coulomb_long * g2 
      VecPointwiseMult(BHD->v[dim], BHD->g2O, BHD->fO_l[dim]);
      // XXX: F_coulomb_long * g2 - F_coulomb_long
      VecAXPY(BHD->v[dim], -1.0, BHD->fO_l[dim]);
      // XXX: FFT(F_coulomb_long * g2 - F_coulomb_long)
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2OOl_fft[dim], BHD->fft_scratch);
      /* HH */
      VecPointwiseMult(BHD->v[dim], BHD->g2H, BHD->fH[dim]);
      //VecAXPY(BHD->v[dim], -1.0, BHD->fH_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HH_fft[dim], BHD->fft_scratch);
      /* Coulomb long */
      VecPointwiseMult(BHD->v[dim], BHD->g2H, BHD->fH_l[dim]);
      VecAXPY(BHD->v[dim], -1.0, BHD->fH_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HHl_fft[dim], BHD->fft_scratch);
      /* HO */
      VecPointwiseMult(BHD->v[dim], BHD->g2HO, BHD->fHO[dim]);
      //VecAXPY(BHD->v[dim], -1.0, BHD->fHO_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HO_fft[dim], BHD->fft_scratch);
      /* Coulomb long */
      VecPointwiseMult(BHD->v[dim], BHD->g2HO, BHD->fHO_l[dim]);
      VecAXPY(BHD->v[dim], -1.0, BHD->fHO_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HOl_fft[dim], BHD->fft_scratch);


    }





/*   VecView(BHD->v[0],PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

}

void RecomputeInitialSoluteData(State *BHD, real damp, real damp_LJ, real zpad)
{
  DA da;
  PetscScalar ***gHini_vec, ***gOini_vec;
  // PetscScalar ***ucH_vec, ***ucO_vec;
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]);
  real r[3], r_s, h[3], interval[2], beta;
  int x[3], n[3], i[3];
  real epsilonH, epsilonO, epsilonHO;
  real sigmaH, sigmaO, sigmaHO;
  real q2H, q2O, q2HO;

  epsilonH = BHD->LJ_paramsH[0];
  epsilonO = BHD->LJ_paramsO[0]; /* FIXME: how comes it is unused? */
  epsilonHO = BHD->LJ_paramsHO[0];

  sigmaH = BHD->LJ_paramsH[1];
  sigmaO = BHD->LJ_paramsO[1];  /* FIXME: how comes it is unused? */
  sigmaHO = BHD->LJ_paramsHO[1];

  q2H = BHD->LJ_paramsH[2];
  q2O = BHD->LJ_paramsO[2];     /* FIXME: how comes it is unused? */
  q2HO = BHD->LJ_paramsHO[2];


  const ProblemData *PD = BHD->PD;
  da = BHD->da;

  PetscPrintf(PETSC_COMM_WORLD,"Recomputing solute data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);


  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));



  VecSet(BHD->g_ini[0], 0.0);
  VecSet(BHD->g_ini[1], 0.0);
  VecSet(BHD->gHO_ini, 0.0);
  VecSet(BHD->uc[0], 0.0);
  VecSet(BHD->uc[1], 0.0);
  VecSet(BHD->ucHO, 0.0);
  FOR_DIM
    {
      VecSet(BHD->fH_l[dim],0.0);
      VecSet(BHD->fO_l[dim],0.0);
      VecSet(BHD->fHO_l[dim],0.0);
    }


  DAVecGetArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecGetArray(da, BHD->g_ini[1], &gOini_vec);

/*   DAVecGetArray(da, BHD->uc[0], &ucH_vec); */
/*   DAVecGetArray(da, BHD->uc[1], &ucO_vec); */
  FOR_DIM
    {
      DAVecGetArray(da, BHD->fH_l[dim], &(fHl_vec[dim]));
      DAVecGetArray(da, BHD->fO_l[dim], &(fOl_vec[dim]));
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

          /*
           * FIXME: it  looks very  strange that gOini_vec  is defined
           * using epsilonHO,  sigmaHO, and q2HO  whereas gHini_vec is
           * defined using epsilonH, sigmaH, and q2H. Is this a bug?
           */

          /* Lennard-Jones */
          gHini_vec[i[2]][i[1]][i[0]] +=
            damp_LJ * beta* Lennard_Jones( r_s, epsilonH, sigmaH);
          gOini_vec[i[2]][i[1]][i[0]] +=
            damp_LJ * beta* Lennard_Jones( r_s, epsilonHO, sigmaHO);

          /* Coulomb short */
          gHini_vec[i[2]][i[1]][i[0]] +=
            damp*beta* Coulomb_short( r_s, q2H);
          gOini_vec[i[2]][i[1]][i[0]] +=
            damp*beta* Coulomb_short( r_s, q2HO);


          /* Coulomb long */
/*        gHini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsH); */
/*        gOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */
/*        gHOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */

           /* Coulomb */
/*        gHini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsH); */
/*        gOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsO); */
/*        gHOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsHO); */

/*        ucH_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */
/*        ucO_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */



        }

  DAVecRestoreArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecRestoreArray(da, BHD->g_ini[1], &gOini_vec);
/*   DAVecRestoreArray(da, BHD->uc[0], &ucH_vec); */
/*   DAVecRestoreArray(da, BHD->uc[1], &ucO_vec); */

  VecCopy( BHD->ucHO, BHD->uc[1]);

/*   ComputeFFTSoluteII(BHD, BHD->uc[0] , BHD->ucHO, BHD->LJ_paramsHO, damp, zpad); */
/*   VecScale(BHD->uc[0], beta); */
/*   VecAXPY( BHD->g_ini[0], beta, BHD->ucHO); */

/*   ComputeFFTSoluteII(BHD, BHD->uc[1] , BHD->ucHO, BHD->LJ_paramsO,  damp, zpad); */
/*   VecAXPY( BHD->g_ini[1], beta, BHD->ucHO); */
/*   VecScale(BHD->uc[1], beta); */

  /* Shift uc's */
/*   VecSum(BHD->uc[0], &fac); */
/*   VecShift(BHD->uc[0], -fac/N[0]/N[1]/N[2]); */
/*   VecSum(BHD->uc[1], &fac); */
/*   VecShift(BHD->uc[1], -fac/N[0]/N[1]/N[2]); */

/*   VecView(BHD->g_ini[0],PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

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
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

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
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

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

typedef struct StepDataStruct
{
  State *BHD ;
  Vec dg_newO, dgO, dgH;
}*StepData;

static PetscErrorCode ComputeStepFunction(SNES snes, Vec x, Vec f, void *data)
{
  StepData SD;
  real con, sumO, sumH; // res
  PetscScalar *x_vec, *f_vec;
  State *BHD;
  Vec gO, gH, dgO2;

  SD = (StepData) data;
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
  StepData SD;
  real step;

  SD = (StepData) malloc(sizeof(*SD));
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
 * This function is the main entry point for the BGY3dM equation for a
 * 2-site solvent and an arbitrary solute.  I guess H2O in the name is
 * a historical baggage.
 */
Vec BGY3dM_solve_H2O_2site(const ProblemData *PD, Vec g_ini)
{
  real norm;
  Vec t_vec;                    /* used for all sites */
  Vec  g[2], dg[2], dg_acc, work;
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

  /* Solutes index */
  const int solute = PD->solute; /* HCl by default */

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

      DACreateGlobalVector (BHD.da, &(g[i]));
      DACreateGlobalVector (BHD.da, &(dg[i]));
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

  /* XXX: Here g0 = beta  * (VM_LJ + VM_coulomb_short) actually.  See:
          (5.106) and (5.108) in Jager's thesis. It is not filled with
          data yet, I assume. */
  const Vec *g0 = BHD.g_ini;    /* just an alias */

  /* set initial guess*/
  VecSet(dg[0],0);
  VecSet(dg[1],0);

  if (bgy3d_getopt_test ("--from-g2")) {
      ComputedgFromg (dg[0], g0[0], BHD.g2HO);
      ComputedgFromg (dg[1], g0[1], BHD.g2O);
  }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      bgy3d_load_vec ("dg0.bin", &dg[0]); /* dgH */
      bgy3d_load_vec ("dg1.bin", &dg[1]); /* dgO */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  }

  for (real damp = damp_start; damp <= 1; damp += 0.1) {

      /* FIXME: I  guess the  logic with damping  factors can  be made
         more straightforward:  */
      real damp_LJ = (damp >= 0 ? 1.0 : 0.0); /* yes, >=, so the
                                                 original */

      /* XXX: Return F * g2.  Note the calculation of F is divided due
              to long  range Coulomb interation.  See  comments in the
              function.  Here F is force within solvents particles.

              In  RecomputeInitialFFTs(),   BHD.uc[0],  BHD.uc[1]  and
              BHD.ucHO  also get  updated  by ComputeFFTfromCoulomb(),
              but  BHD.uc[0]   and  BHD.uc[1]  are   re-calculated  by
              RecomputeInitialSoluteData_XXX  again  and BHD.ucHO  are
              not used

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

              and sum to the solution by the end of each iteration for
              solvent site  A and B,  in the code hereafter,  they are
              stored as BHD.uc[0] and BHD.uc[1] */

      RecomputeInitialFFTs(&BHD, (damp > 0.0 ? damp : 0.0), 1.0);

      /* FIXME:  avoid  storing  vectors  fg2XY* on  the  grid  across
         iterations. Only a scalr kernel is needed: */
      kernel (BHD.da, BHD.PD, BHD.fg2HH_fft, NULL, ker_fft_S[0][0]);
      kernel (BHD.da, BHD.PD, BHD.fg2HO_fft, NULL, ker_fft_S[0][1]); /* == [1][0] */
      kernel (BHD.da, BHD.PD, BHD.fg2OO_fft, NULL, ker_fft_S[1][1]);

      kernel (BHD.da, BHD.PD, BHD.fg2HHl_fft, BHD.ucH_fft,  ker_fft_L[0][0]);
      kernel (BHD.da, BHD.PD, BHD.fg2HOl_fft, BHD.ucHO_fft, ker_fft_L[0][1]); /* == [1][0] */
      kernel (BHD.da, BHD.PD, BHD.fg2OOl_fft, BHD.ucO_fft,  ker_fft_L[1][1]);

      /* FIXME: what is  the point to split the  kernel in two pieces?
         Redefine S := S + L and forget about L */
      for (int i = 0; i < 2; i++)
          for (int j = 0; j <= i; j++) {
              /* S := (damp / damp0) * L + damp_LJ * S */
              bgy3d_fft_axpby (BHD.da, ker_fft_S[i][j], damp / damp0, damp_LJ, ker_fft_L[i][j]);
          }

      /* XXX: Return BHD.g_ini[0], BHD.g_ini[1] (see definition above)
              and BHD.uc[0], BHD.uc[1], which are VM_Coulomb_long, but
              should they  multiply by  beta?  Solute is  hardcoded as
              HCl (==0) for standard test */
      bgy3d_solute_field (&BHD, solute, (damp > 0.0 ? damp : 0.0), 1.0);

      /* FIXME:  Check if this  is redundant  --- it  was mechanically
         moved from  the body of the  above func (because  it does not
         belong into solute code): */
      FOR_DIM {
        VecSet(BHD.fH_l[dim], 0.0);
        VecSet(BHD.fO_l[dim], 0.0);
        VecSet(BHD.fHO_l[dim], 0.0); /* What is it used for? */
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
             RecomputeInitialFFTs(). FIXME: avoid this ugliness.

             I assume  BHD.uc[] was not  modified since it was  set in
             the  solute   specific  RecomputeInitialSoluteData()  and
             before it first use  in VecAXPY() below.  Note that uc[0]
             and  uc[1] differ from  the true  Coulomb potential  by a
             factor equal  to the charge of the  respective site.  Why
             do  we   keep  two  versions  of   essentially  the  same
             potential?  Why not adapt the factor here instead? */

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

              /* Add Coulomb field BHD.uc[i] to: */
              VecAXPY(dg_acc, 1.0, BHD.uc[i]);

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

          /* Now  that dg[]  has bee  computed one  can  safely update
             g[]: */
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

  /* Save final distribution, use binary format: */
  bgy3d_save_vec ("g0.bin", g[0]); /* gH */
  bgy3d_save_vec ("g1.bin", g[1]); /* gO */

  /* Clean up and exit ... */
  bgy3d_fft_free (dg_acc_fft);

  for (int i = 0; i < 2; i++) {
      VecDestroy (g[i]);
      VecDestroy (dg[i]);

      bgy3d_fft_free (g_fft[i]);

      for (int j = 0; j <= i; j++) {
          bgy3d_fft_free (ker_fft_S[i][j]);
          bgy3d_fft_free (ker_fft_L[i][j]);
      }
  }

  VecDestroy(dg_acc);
  VecDestroy(work);
  VecDestroy(t_vec);

  finalize_state (&BHD);

  return PETSC_NULL;
}


Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini)
{
  real a1, a, damp, damp_LJ;
  real count=0.0, norm, aO;
  int iter;
  Vec g0H, g0O, dgH, dgO,  dg_new, dg_new2, f, gH, gO;
  Vec tH, tO, dg_newH, dg_newO;
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

  DACreateGlobalVector(BHD.da, &dg_histO);
  DACreateGlobalVector(BHD.da, &dg_histH);
  VecSet(dg_histH, 0.0);
  VecSet(dg_histO, 0.0);

  g0H=BHD.g_ini[0];
  g0O=BHD.g_ini[1];

/*   VecView(BHD.fHO_l[0],PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dg_new,0.0);

  if (bgy3d_getopt_test ("--from-g2")) {
      ComputedgFromg (dgH, g0H, BHD.g2HO);
      ComputedgFromg (dgO, g0O, BHD.g2O);
  }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      bgy3d_load_vec ("dg0.bin", &dgH); /* dgH */
      bgy3d_load_vec ("dg1.bin", &dgO); /* dgO */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  }

  for( damp=damp_start; damp <=1; damp+=0.1)
    {
      if(damp==-0.01)
        {
          damp_LJ=0;
          //a0=0.4;
          RecomputeInitialFFTs(&BHD, 0.0, 1.0);
          bgy3d_solute_field (&BHD, /* Butanoic Acid */ 4, 0.0, 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }
      else if(damp==0.0)
        {
          damp_LJ=1.0;
          //a0=0.5;
          RecomputeInitialFFTs(&BHD, 0.0, 1.0);
          bgy3d_solute_field (&BHD, /* Butanoic Acid */ 4, 0.0, 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }
      else
        {
          damp_LJ=1.0;
          //a0=0.0002/damp;
          count+=1.0;
          //a0 = 0.1/4./(count);
          a0 = 0.1/(count+5.0);
          RecomputeInitialFFTs(&BHD, (damp), 1.0);
          bgy3d_solute_field (&BHD, /* Butanoic Acid */ 4, damp, 1.0);
          PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
        }

      /* FIXME:  Check if this  is redundant  --- it  was mechanically
         moved from the body of the bgy3d_solute_field() func (because
         it does not belong into solute code): */
      FOR_DIM {
        VecSet(BHD.fH_l[dim], 0.0);
        VecSet(BHD.fO_l[dim], 0.0);
        VecSet(BHD.fHO_l[dim], 0.0); /* What is it used for? */
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


  /*     VecCopy(BHD.g2HO, gH);  */
/*       VecCopy(BHD.g2O, gO);  */

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
          Compute_H2O_interS(&BHD, BHD.fg2HO_fft, gO, BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.fg2HH_fft, gH, BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fg2HOl_fft, gO, BHD.ucHO_fft, (damp / damp0) * BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fg2HHl_fft, gH, BHD.ucH_fft, (damp / damp0) * BHD.rhos[0], dg_new2);
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


          //VecShift(dg_new, -0.01);
          //VecAXPY(dg_new, dgH_norm,  BHD.uc[0]);
          VecAXPY(dg_new, 1.0, BHD.uc[0]);
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
          Compute_H2O_interS(&BHD, BHD.fg2OO_fft, gO, BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.fg2HO_fft, gH, BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fg2OOl_fft, gO, BHD.ucO_fft, (damp / damp0) * BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fg2HOl_fft, gH, BHD.ucHO_fft, (damp / damp0) * BHD.rhos[0], dg_new2);
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


          //VecShift(dg_new, 0.01);
          //VecAXPY(dg_new, dgO_norm, BHD.uc[1]);
          VecAXPY(dg_new, 1, BHD.uc[1]);
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
  VecDestroy(dg_newH);
  VecDestroy(dg_newO);
  VecDestroy(dg_histH);
  VecDestroy(dg_histO);

  finalize_state (&BHD);


  //VecView(BHD.uc[1],PETSC_VIEWER_STDERR_WORLD);
   //exit(1);

  return PETSC_NULL;
}
