/*==========================================================*/
/*  $Id: bgy3dH2OS.c,v 1.20 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"

#include "bgy3d_SolventParameters.h"
#include "bgy3d-getopt.h"
#include "bgy3d-solutes.h"
#include "bgy3dH2OS.h"

#define damp0 0.001

extern real NORM_REG;

static State *BGY3dH2OData_malloc(ProblemData *PD)
{
  State *BHD;
  DA da;
  real beta, maxL;
  int x[3], n[3];
  int np;
  int local_nx, local_x_start, local_ny, local_y_start, total_local_size;
  PetscInt lx[1], ly[1], *lz;
  PetscErrorCode ierr;

  BHD = (BGY3dH2OData) malloc(sizeof(*BHD));

  /****************************************************/
  /* set Lennard-Jones and Coulomb parameters */
  /****************************************************/

  /* water hydrogen */
  BHD->LJ_paramsH[0] = eH;  /* epsilon  */
  BHD->LJ_paramsH[1] = sH;  /* sigma    */
  BHD->LJ_paramsH[2] = SQR(qH); /* charge product */

  /* water oxygen */
  BHD->LJ_paramsO[0] = eO;  /* epsilon  */
  BHD->LJ_paramsO[1] = sO;  /* sigma    */
  BHD->LJ_paramsO[2] = SQR(qO); /* charge product */

  /* water O-H mixed parameters */
  BHD->LJ_paramsHO[0] = sqrt(eH*eO);  /* epsilon  */
  BHD->LJ_paramsHO[1] = 0.5*(sH+sO);  /* sigma    */
  BHD->LJ_paramsHO[2] = qH*qO; /* charge product */

  /****************************************************/



  BHD->PD = PD;
  /*****************************/
  /* reset standard parameters */
  /*****************************/
  maxL=12.0;
  bgy3d_getopt_real ("-L", &maxL);
  PD->interval[0] = -maxL;//-25.0;
  PD->interval[1] = maxL;//25.0;
  FOR_DIM
    PD->h[dim] = (PD->interval[1]-PD->interval[0])/PD->N[dim];
  PD->N3 = PD->N[0]*PD->N[1]*PD->N[2];
  beta = 1.6889;
  bgy3d_getopt_real ("-beta", &beta);
  PD->beta = beta;
  PetscPrintf(PETSC_COMM_WORLD, "Corrected domain size:\n");
  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  //PetscPrintf(PETSC_COMM_WORLD, "Boundary smoothing parameters : SL= %f  SR= %f\n", SL, SR);
  //PetscPrintf(PETSC_COMM_WORLD, "ZEROPAD= %f\n", ZEROPAD);
  //PetscPrintf(PETSC_COMM_WORLD, "Regularization of normalization: NORM_REG= %e\n", NORM_REG);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", PD->beta);
  /******************************/
  BHD->beta = PD->beta;
  BHD->rho  = PD->rho;
  BHD->rho_H = PD->rho;
  BHD->rho_O = PD->rho;

  /* Initialize parallel stuff: fftw + petsc */
  BHD->fft_plan_fw = fftw3d_mpi_create_plan(PETSC_COMM_WORLD,
					    PD->N[2], PD->N[1], PD->N[0],
					    FFTW_FORWARD, FFTW_ESTIMATE);
  BHD->fft_plan_bw = fftw3d_mpi_create_plan(PETSC_COMM_WORLD,
					    PD->N[2], PD->N[1], PD->N[0],
					    FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwnd_mpi_local_sizes(BHD->fft_plan_fw, &local_nx, &local_x_start,
			 &local_ny, &local_y_start, &total_local_size);
  /* Get number of processes */
  MPI_Comm_size(PETSC_COMM_WORLD, &np);

  /* Create Petsc Distributed Array according to fftw data distribution*/
  lz = (PetscInt*) malloc(np*sizeof(*lz));

  MPI_Allgather( &local_nx, 1, MPI_INT, lz, 1, MPI_INT, PETSC_COMM_WORLD);
  ly[0]=PD->N[1];
  lx[0]=PD->N[2];

#ifdef L_BOUNDARY
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     1, 1, np,
	     1,1,
	     lx, ly, lz,
	     &(BHD->da));
  da = BHD->da;
   /* Create Matrix with appropriate non-zero structure */
  DAGetMatrix( da, MATMPIAIJ, &(BHD->M));
  DACreateGlobalVector(da, &(BHD->xH));
  DACreateGlobalVector(da, &(BHD->xO));
  VecSet(BHD->xH, 0.0);
  VecSet(BHD->xO, 0.0);
#endif
#ifdef L_BOUNDARY_MG
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     1, 1, np,
	     1,1,
	     lx, ly, lz,
	     &(BHD->da));
  da = BHD->da;

  for(dim=0; dim<np; dim++)
    lz[dim]/=2;
  lx[0]/=2;
  ly[0]/=2;
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0]/2, PD->N[1]/2, PD->N[2]/2,
	     1, 1, np,
	     1,1,
	     lx, ly, lz,
	     &(BHD->da_dmmg));
#endif

#ifndef  L_BOUNDARY
#ifndef  L_BOUNDARY_MG
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     1, 1, np,
	     1,0,
	     lx, ly, lz,
	     &(BHD->da));


  da = BHD->da;
#endif
#endif

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

  /*
   * CHKERRQ()  is  a  macro that  expands  to  an  if with  a  return
   * statement  in  the if-block.   It  is  not  suitable for  use  in
   * functions  returning  anything   that  but  PetscErrorCode.  This
   * function returns BGY3dH2OData.
   */

  /* Create global vectors */
  ierr = DACreateGlobalVector(da, &(BHD->g_ini[0])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g_ini[1])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->gHO_ini)); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->uc[0])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->uc[1])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->ucHO)); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g2H)); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g2O)); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g2HO)); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->pre)); // CHKERRQ(ierr);
  assert (!ierr);

  FOR_DIM
    {

      ierr = DACreateGlobalVector(da, &(BHD->fH[dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->fO[dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->fHO[dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->fH_l[dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->fO_l[dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->fHO_l[dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->v[dim])); // CHKERRQ(ierr);
      assert (!ierr);
    }


  VecSet(BHD->pre,0.0);

/*   VecView(BHD->gHO_ini,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */



  if(BHD->fft_plan_fw == NULL || BHD->fft_plan_bw == NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Failed to get fft_plan of proc %d.\n",
		  PD->id);
      exit(1);
    }


  /* Allocate memory for fft */
  FOR_DIM
    {
      BHD->fg2HH_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
      BHD->fg2OO_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
      BHD->fg2HO_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
      BHD->fg2HHl_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
      BHD->fg2OOl_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
      BHD->fg2HOl_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
      BHD->fg2_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));

      BHD->fO_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
      BHD->fH_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
    }

  BHD->g_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->gfg2_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->fft_scratch = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->ucH_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->ucO_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->ucHO_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->wHO_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->wHH_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));



  /* Read g^2  from file */
/*   ReadPairDistribution(BHD, "g2_OO", BHD->g2O); */
/*   ReadPairDistribution(BHD, "g2_HH", BHD->g2H); */
/*   ReadPairDistribution(BHD, "g2_HO", BHD->g2HO); */

#ifdef CS2
  ReadPairDistribution(BHD, "g2C", BHD->g2O);
  ReadPairDistribution(BHD, "g2S", BHD->g2H);
  ReadPairDistribution(BHD, "g2CS", BHD->g2HO);
#else
  bgy3d_load_vec ("g2H.bin", &(BHD->g2H));
  bgy3d_load_vec ("g2O.bin", &(BHD->g2O));
  bgy3d_load_vec ("g2HO.bin", &(BHD->g2HO));
#endif






  /* Compute initial data */
  //RecomputeInitialFFTs(BHD, 1.0, 1.0);

  /* Compute Solute dependent initial data */
  //RecomputeInitialSoluteData(BHD, 1.0, 1.0);

  free(lz);

  return BHD;
}



static void BGY3dH2OData_free2(State *BHD)
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
      free(BHD->fg2HH_fft[dim]);
      free(BHD->fg2OO_fft[dim]);
      free(BHD->fg2HO_fft[dim]);
      free(BHD->fg2HHl_fft[dim]);
      free(BHD->fg2OOl_fft[dim]);
      free(BHD->fg2HOl_fft[dim]);
      free(BHD->fg2_fft[dim]);

      free(BHD->fO_fft[dim]);
      free(BHD->fH_fft[dim]);
    }
  free(BHD->g_fft);
  free(BHD->gfg2_fft);
  free(BHD->fft_scratch);
  free(BHD->ucH_fft);
  free(BHD->ucO_fft);
  free(BHD->ucHO_fft);
  free(BHD->wHO_fft);
  free(BHD->wHH_fft);

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
  VecDestroy(BHD->xO);
  VecDestroy(BHD->xH);
#endif
#ifdef L_BOUNDARY_MG
  DMMGDestroy(BHD->dmmg);
#endif
  DADestroy(BHD->da);

  fftwnd_mpi_destroy_plan(BHD->fft_plan_fw);
  fftwnd_mpi_destroy_plan(BHD->fft_plan_bw);

  free(BHD);
}

#ifdef L_BOUNDARY
/* Initialize M-Matrix with appropriate stencil */
void InitializeLaplaceMatrix(State *BHD, real zpad)
{
  ProblemData *PD;
  Mat M;
  DA da;
  int x[3], n[3], i[3], N[3], border;
  MatStencil col[3],row;
  PetscScalar v[3], vb=1.0;
  real h[3];

  PetscPrintf(PETSC_COMM_WORLD,"Assembling Matrix...");

  da = BHD->da;
  PD = BHD->PD;
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

static void CopyBoundary(State *BHD, Vec gfrom, Vec gto, real zpad)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], N[3], border; // ic[3], k;
  PetscScalar ***gfrom_vec, ***gto_vec;

  PD = BHD->PD;
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
/* 	for( i[ic[1]]=x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++) */
/* 	  for( i[ic[2]]=x[ic[2]]; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++) */
/* 	    gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]]; */
/*       i[ic[0]] = N[ic[0]]-1; */
/*       if( x[ic[0]]+n[ic[0]] == N[ic[0]] )  */
/* 	for( i[ic[1]]=x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++) */
/* 	  for( i[ic[2]]=x[ic[2]]; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++) */
/* 	    gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]]; */
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

real ImposeLaplaceBoundary(State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter)
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

void ReadPairDistribution(State *BHD, char *filename, Vec g2)
{
  ProblemData *PD;
  DA da;
  FILE *fp;
  real *xg, *g;
  real r[3], r_s, h[3], interval[2];
  int index=0;
  int x[3], n[3], i[3], k;
  PetscScalar ***g2_vec;

  da = BHD->da;
  PD = BHD->PD;
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
  ProblemData *PD;
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


  PD = BHD->PD;
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
/* 	  gH_vec[i[2]][i[1]][i[0]] +=  */
/* 	    exp(-damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsH)); */
/* 	  gO_vec[i[2]][i[1]][i[0]] +=  */
/* 	    exp(-damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsO)); */
/* 	  gHO_vec[i[2]][i[1]][i[0]] +=  */
/* 	    exp(-damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsHO)); */


	  /* omega_HH + omega_HO */
/* 	  wHO_vec[i[2]][i[1]][i[0]] = wconst_HO * exp(-SQR(wG)*SQR(r_s-r_HO)); */
/* 	  wHH_vec[i[2]][i[1]][i[0]] = wconst_HH * exp(-SQR(wG)*SQR(r_s-r_HH)); */
/* 	  wHO_vec[i[2]][i[1]][i[0]] = damp * Coulomb_long( r_s, BHD->LJ_paramsO); */

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
/*  	      fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

	      /* Coulomb long */
 /*  	      fHl_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fOl_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHOl_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

	      /* Coulomb */
/* 	      fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsHO); */


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
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2OO_fft[dim],
			     BHD->fft_scratch, x, n, 0);
      /* Coulomb long */
      // XXX: F_coulomb_long * g2 
      VecPointwiseMult(BHD->v[dim], BHD->g2O, BHD->fO_l[dim]);
      // XXX: F_coulomb_long * g2 - F_coulomb_long
      VecAXPY(BHD->v[dim], -1.0, BHD->fO_l[dim]);
      // XXX: FFT(F_coulomb_long * g2 - F_coulomb_long)
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2OOl_fft[dim],
			     BHD->fft_scratch, x, n, 0);
      /* HH */
      VecPointwiseMult(BHD->v[dim], BHD->g2H, BHD->fH[dim]);
      //VecAXPY(BHD->v[dim], -1.0, BHD->fH_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HH_fft[dim],
			     BHD->fft_scratch, x, n, 0);
      /* Coulomb long */
      VecPointwiseMult(BHD->v[dim], BHD->g2H, BHD->fH_l[dim]);
      VecAXPY(BHD->v[dim], -1.0, BHD->fH_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HHl_fft[dim],
			     BHD->fft_scratch, x, n, 0);
      /* HO */
      VecPointwiseMult(BHD->v[dim], BHD->g2HO, BHD->fHO[dim]);
      //VecAXPY(BHD->v[dim], -1.0, BHD->fHO_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HO_fft[dim],
			     BHD->fft_scratch, x, n, 0);
      /* Coulomb long */
      VecPointwiseMult(BHD->v[dim], BHD->g2HO, BHD->fHO_l[dim]);
      VecAXPY(BHD->v[dim], -1.0, BHD->fHO_l[dim]);
      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], BHD->fg2HOl_fft[dim],
			     BHD->fft_scratch, x, n, 0);


    }





/*   VecView(BHD->v[0],PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

}

void RecomputeInitialSoluteData(State *BHD, real damp, real damp_LJ, real zpad)
{
  DA da;
  ProblemData *PD;
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


  PD = BHD->PD;
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
/* 	  gHini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsH); */
/* 	  gOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */
/* 	  gHOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */

	   /* Coulomb */
/* 	  gHini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb( r_s, BHD->LJ_paramsH); */
/* 	  gOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb( r_s, BHD->LJ_paramsO); */
/* 	  gHOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb( r_s, BHD->LJ_paramsHO); */

/* 	  ucH_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */
/* 	  ucO_vec[i[2]][i[1]][i[0]] +=  */
/* 	    damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */



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


void Compute_H2O_interS(State *BHD,
			fftw_complex *(fg2_fft[3]), Vec g, fftw_complex *coul_fft,
			fftw_complex *(fs_fft[3]), real con, real rho, Vec dg_help)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  fftw_complex *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, sign; // confac;

  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];


  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  fac = L/(2.*M_PI);  /* BHD->f ist nur grad U, nicht F=-grad U  */

  /* confac = SQR(M_PI/L/2.); */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* rho*F*g^2 g*/
  /************************************************/


  /* fft(g) */

  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g, g_fft, scratch,
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
	      dg_fft[index].re = 0;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = h*h*fac/k;
	      /* phase shift factor for x=x+L/2 */
	      sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);


	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac * sign *
		(fg2_fft[dim][index].re * g_fft[index].im
		 + fg2_fft[dim][index].im * g_fft[index].re) ;


	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac * sign *
		(-fg2_fft[dim][index].re * g_fft[index].re
		 + fg2_fft[dim][index].im * g_fft[index].im);




	      /* right one: */
	      /*****************************/
	      /* long range Coulomb part */
/* 	      dg_fft[index].re += h*sign* */
/* 		(coul_fft[index].re*(g_fft[index].re + 0*con*sign/(h*rho)) */
/* 		 - coul_fft[index].im*g_fft[index].im); */

/* 	      dg_fft[index].im += h*sign* */
/* 		(coul_fft[index].re*g_fft[index].im */
/* 		 + coul_fft[index].im*(g_fft[index].re + 0*con*sign/(h*rho))); */
	      /******************************/

	      /*****************************/
	      /* long range Coulomb part */
/* 	      FOR_DIM */
/* 		dg_fft[index].re += ic[dim] * k_fac *  */
/* 		( fs_fft[dim][index].im * con/(h*rho)) ; */

/* 	      FOR_DIM  */
/* 		dg_fft[index].im += ic[dim] * k_fac *  */
/* 		(-fs_fft[dim][index].re * con/(h*rho)); */
	      /******************************/


	      /*****************************/
	      /* long range Coulomb part */
/* 	      dg_fft[index].re += h*sign* */
/* 		(coul_fft[index].re*(g_fft[index].re + con*sign/(h*rho)*exp(-k*confac)) */
/* 		- coul_fft[index].im*g_fft[index].im); */

/* 	      dg_fft[index].im += h*sign* */
/* 		(coul_fft[index].re*g_fft[index].im */
/* 		 + coul_fft[index].im*(g_fft[index].re + con*sign/(h*rho)*exp(-k*confac))); */
	      /******************************/

	      //dg_fft[index].re = sign*exp(-k*confac);

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */
	    }
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft,
			 scratch, x, n, 0.0);

  VecScale(dg_help, rho*PD->beta/L/L/L);



/*   ImposeLaplaceBoundary(BHD, dg_help, BHD->v[0], BHD->v[1], 8.0); */
/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);   */
/*    exit(1);      */

}

static void Compute_H2O_interS_C(State *BHD,
			  fftw_complex *(fg2_fft[3]), Vec g, fftw_complex *coul_fft,
			  fftw_complex *(fs_fft[3]), real con, real rho, Vec dg_help, real damp)
{
  ProblemData *PD;
  DA da;
  int x[3], n[3], i[3], index, N[3], ic[3];
  fftw_complex *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, sign, dampfac; // confac;

  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];


  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  fac = L/(2.*M_PI);  /* BHD->f ist nur grad U, nicht F=-grad U  */

  /* confac = SQR(M_PI/L/2.); */

  dampfac = damp/damp0;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* rho*F*g^2 g*/
  /************************************************/


  /* fft(g) */

  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g, g_fft, scratch,
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
	      dg_fft[index].re = 0;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = h*h*fac/k;
	      /* phase shift factor for x=x+L/2 */
	      sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);


	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac * sign *
		(fg2_fft[dim][index].re * g_fft[index].im
		 + fg2_fft[dim][index].im * g_fft[index].re) ;


	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac * sign *
		(-fg2_fft[dim][index].re * g_fft[index].re
		 + fg2_fft[dim][index].im * g_fft[index].im);




	      /* right one: */
	      /*****************************/
	      /* long range Coulomb part */
	      dg_fft[index].re += h*sign*
		(coul_fft[index].re*g_fft[index].re
		 - coul_fft[index].im*g_fft[index].im);

	      dg_fft[index].im += h*sign*
		(coul_fft[index].re*g_fft[index].im
		 + coul_fft[index].im*g_fft[index].re );
	      /******************************/

	      /*****************************/
	      /* long range Coulomb part */
/* 	      FOR_DIM */
/* 		dg_fft[index].re += ic[dim] * k_fac *  */
/* 		( fs_fft[dim][index].im * con/(h*rho)) ; */

/* 	      FOR_DIM  */
/* 		dg_fft[index].im += ic[dim] * k_fac *  */
/* 		(-fs_fft[dim][index].re * con/(h*rho)); */
	      /******************************/


	      /*****************************/
	      /* long range Coulomb part */
/* 	      dg_fft[index].re += h*sign* */
/* 		(coul_fft[index].re*(g_fft[index].re + con*sign/(h*rho)*exp(-k*confac)) */
/* 		- coul_fft[index].im*g_fft[index].im); */

/* 	      dg_fft[index].im += h*sign* */
/* 		(coul_fft[index].re*g_fft[index].im */
/* 		 + coul_fft[index].im*(g_fft[index].re + con*sign/(h*rho)*exp(-k*confac))); */
	      /******************************/

	      //dg_fft[index].re = sign*exp(-k*confac);

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */
	    }
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft,
			 scratch, x, n, 0.0);

  VecScale(dg_help, dampfac*rho*PD->beta/L/L/L);



/*   ImposeLaplaceBoundary(BHD, dg_help, BHD->v[0], BHD->v[1], 8.0); */
/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);   */
/*    exit(1);      */

}



real ComputeCharge(State *BHD, Vec g1, Vec g2)
{
  ProblemData *PD;
  real g1_sum, g2_sum, c;
  Vec help;


  PD = BHD->PD;
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
  ProblemData *PD;
  Vec gO, gH, dgO2;

  SD = (StepData) data;
  BHD = SD->BHD;
  PD = BHD->PD;
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
 * This function is the main entry point for the BGY3dM equation for a
 * 2-site solvent and an arbitrary solute.  I guess H2O in the name is
 * a historical baggage.
 */
Vec BGY3dM_solve_H2O_2site(ProblemData *PD, Vec g_ini, int vdim)
{
  State *BHD;
  real a0=0.1, a1, a, damp_start=0.0, norm_tol=1.0e-2, zpad=1000.0, damp, damp_LJ;
  real count=0.0, norm, aO;
  int max_iter=100, iter;
  Vec g0H, g0O, dgH, dgO,  dg_new, dg_new2, f, gH, gO;
  Vec tH, tO, dg_newH, dg_newO;
  PetscScalar dgH_norm, dgO_norm;
  real dgH_old, dgO_old;
  int mycount=0, upwards, namecount=0;
  char nameH[20], nameO[20];



  Vec dg_histO, dg_histH;

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (2-site) equation ...\n");

  BHD = BGY3dH2OData_malloc(PD);
  if( r_HH >0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"WARNING: Solvent not a 2-Site model!\n");
    }

  /* read BGY3d specific things from command line */
  /* Mixing parameter */
  bgy3d_getopt_real ("-lambda", &a0);
  if(a0>1 || a0<0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lambda out of range: lambda=%f\n",a0);
      exit(1);
    }
   /* Get damp_start from command line*/
  bgy3d_getopt_real ("-damp_start", &damp_start);
  /* Number of total iterations */
  bgy3d_getopt_int ("-max_iter", &max_iter);
  /* norm_tol for convergence test */
  bgy3d_getopt_real ("-norm_tol", &norm_tol);
  /* Zeropad */
  bgy3d_getopt_real ("-zpad", &zpad);
  /*********************************/

  PetscPrintf(PETSC_COMM_WORLD,"lambda = %f\n",a0);
  PetscPrintf(PETSC_COMM_WORLD,"tolerance = %e\n",norm_tol);
  PetscPrintf(PETSC_COMM_WORLD,"zpad = %f\n",zpad);
  PetscPrintf(PETSC_COMM_WORLD,"max_iter = %d\n",max_iter);

  ImposeBoundaryCondition_Initialize( BHD, zpad);
#ifdef L_BOUNDARY
  BHD->zpad=zpad;
  /* Assemble Laplacian matrix */
  InitializeLaplaceMatrix(BHD, zpad);
  /* Create KSP environment */
  InitializeKSPSolver(BHD);
#endif
#ifdef L_BOUNDARY_MG
  BHD->zpad=zpad;

  InitializeDMMGSolver(BHD);
#endif
  DACreateGlobalVector(BHD->da, &gH);
  DACreateGlobalVector(BHD->da, &gO);
  DACreateGlobalVector(BHD->da, &dgH);
  DACreateGlobalVector(BHD->da, &dgO);
  DACreateGlobalVector(BHD->da, &dg_new);
  DACreateGlobalVector(BHD->da, &dg_new2);
  DACreateGlobalVector(BHD->da, &f);

  DACreateGlobalVector(BHD->da, &tH);
  DACreateGlobalVector(BHD->da, &tO);

  DACreateGlobalVector(BHD->da, &dg_newH);
  DACreateGlobalVector(BHD->da, &dg_newO);

  DACreateGlobalVector(BHD->da, &dg_histO);
  DACreateGlobalVector(BHD->da, &dg_histH);
  VecSet(dg_histH, 0.0);
  VecSet(dg_histO, 0.0);

// XXX: here g0 = beta * (VM_LJ + VM_coulomb_short) actually
// XXX: see: (5.106) and (5.108) in Jager's thesis
  g0H=BHD->g_ini[0];
  g0O=BHD->g_ini[1];

/*   VecView(BHD->fHO_l[0],PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dg_new,0.0);

  if (bgy3d_getopt_test ("-fromg2")) {
      ComputedgFromg (dgH, g0H, BHD->g2HO);
      ComputedgFromg (dgO, g0O, BHD->g2O);
  }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("-loadH2O")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      bgy3d_load_vec ("dgH.bin", &dgH); /* dgH */
      bgy3d_load_vec ("dgO.bin", &dgO); /* dgO */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  }

  for (damp = damp_start; damp <= 1; damp += 0.1) {

      /* FIXME: I  guess the  logic with damping  factors can  be made
         more straightforward:  */
      if (damp > 0)
	  count += 1.0;

      damp_LJ = (damp >= 0 ? 1.0 : 0.0); /* yes, >=, so the
                                            original */

      /* XXX: Return F * g2.  Note the calculation of F is divided due
              to long  range Coulomb interation.  See  comments in the
              function.  Here F is force within solvents particles. */
      RecomputeInitialFFTs(BHD, (damp > 0.0 ? damp : 0.0), 1.0);

      /* XXX: Return  BHD->g_ini[0],   BHD->g_ini[1]  (see  definition
              above)    and   BHD->uc[0],   BHD->uc[1],    which   are
              VM_Coulomb_long,  but  should  they  multiply  by  beta?
              Solute is hardcoded as HCl for standard test */
      RecomputeInitialSoluteData_HCl(BHD, (damp > 0.0 ? damp : 0.0), 1.0, zpad);

      PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);

      //Smooth_Function(BHD, g0O, zpad-1, zpad, 0.0);
      //Smooth_Function(BHD, g0H, zpad-1, zpad, 0.0);

      /* XXX: See  p116-177 in thesis:  Boundary Conditions  (5.107) -
             (5.110):  first  impose  boundary condistion  then  solve
             laplacian equation and substrate from g0. */
      ImposeLaplaceBoundary(BHD, g0H, tH, BHD->xH, zpad, NULL);
      ImposeLaplaceBoundary(BHD, g0O, tH, BHD->xO, zpad, NULL);

      /* XXX: then correct g0 with boundary condition again */
      Zeropad_Function(BHD, g0O, zpad, 0.0);
      Zeropad_Function(BHD, g0H, zpad, 0.0);
      /* g=g0*exp(-dg) */

      ComputeH2O_g( gH, g0H, dgH);
      ComputeH2O_g( gO, g0O, dgO);


  /*     VecCopy(BHD->g2HO, gH);  */
/*       VecCopy(BHD->g2O, gO);  */

      a=a0;
      a1=a0;
      for(iter=0; iter<max_iter; iter++)
	{

	  if( !(iter%10) && iter>0 )
	    a=a1;
	  else
	    a=a0;

/* 	  if( !(iter%50) && iter>0 && NORM_REG>=5.0e-2) */
/* 	    NORM_REG/=2.; */

	  PetscPrintf(PETSC_COMM_WORLD,"iter %d: function norms: %e ", iter+1, NORM_REG);


	  /* H */
	  VecSet(dg_new,0.0);
	  Compute_H2O_interS(BHD,
			     BHD->fg2HO_fft, gO, BHD->ucHO_fft, BHD->fH_fft,
			     1.0, BHD->rho_O, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS(BHD,
			     BHD->fg2HH_fft, gH, BHD->ucH_fft, BHD->fH_fft,
			     0.0, BHD->rho_H, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  VecScale(dg_new,damp_LJ);

	  /* Coulomb long */
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2HOl_fft, gO, BHD->ucHO_fft, BHD->fH_fft,
			       1.0, BHD->rho_O, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2HHl_fft, gH, BHD->ucH_fft, BHD->fH_fft,
			       0.0, BHD->rho_H, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);


	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gO, tO , dg_new2, f, zpad);
	  //Compute_dg_H2O_intra_lnIII(BHD, gO, tO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gO, tO, r_HO, dg_new2, f);
	  Compute_dg_H2O_intra_ln(BHD, tO, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);

/*  	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);      */
/*  	  exit(1);    */

	  //dgH_norm = ComputeCharge(BHD, gH, gO);
	  //PetscPrintf(PETSC_COMM_WORLD, " %e ", dgH_norm);
/* 	  if(dgH_norm<1.0) */
/* 	    VecScale(dg_new, 1.2); */
/* 	  else if(dgH_norm >1.0) */
/* 	    VecScale(dg_new, 0.8); */


	  //VecShift(dg_new, -0.01);
	  //VecAXPY(dg_new, dgH_norm,  BHD->uc[0]);
	  VecAXPY(dg_new, 1.0, BHD->uc[0]);
/* 	  dgH_norm = ImposeBoundaryConditionII(BHD, dg_new, zpad); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " %e ", dgH_norm); */
	  //VecShift(dg_new, -dgH_norm);
	  //ImposeBoundaryCondition( BHD, dg_new);

	  ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xH, zpad, NULL);
	  Zeropad_Function(BHD, dg_new, zpad, 0.0);
	  //Smooth_Function(BHD, dg_new, zpad-1, zpad, 0.0);




	  VecCopy(dg_new, dg_newH);

/* 	  VecCopy(dgH, f); */
/* 	  VecAXPBY(dgH, a, (1-a), dg_new); */
/* 	  VecAXPY(f, -1.0, dgH); */
/* 	  VecNorm(f, NORM_INFINITY, &dgH_norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"H= %e  ", dgH_norm/a); */
	  //ComputeH2O_g( gH, g0H, dgH);


	  /* O */
	  VecSet(dg_new,0.0);
	  Compute_H2O_interS(BHD,
			     BHD->fg2OO_fft, gO, BHD->ucO_fft, BHD->fO_fft,
			     1.0, BHD->rho_O, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS(BHD,
			     BHD->fg2HO_fft, gH, BHD->ucHO_fft, BHD->fO_fft,
			     0.0, BHD->rho_H, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  VecScale(dg_new,damp_LJ);

	  /* Coulomb long */
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2OOl_fft, gO, BHD->ucO_fft, BHD->fO_fft,
			       1.0, BHD->rho_O, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2HOl_fft, gH, BHD->ucHO_fft, BHD->fO_fft,
			       0.0, BHD->rho_H, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);



	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gH, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gH, tH, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);

/*  	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);      */
/*  	  exit(1);    */


	  //dgO_norm = ComputeCharge(BHD, gH, gO);
	  //PetscPrintf(PETSC_COMM_WORLD, " %e ", dgO_norm);
/* 	  if(dgO_norm<1.0) */
/* 	    VecScale(dg_new, 0.9); */
/* 	  else if(dgO_norm >1.0) */
/* 	    VecScale(dg_new, 1.1); */


	  //VecShift(dg_new, 0.01);
	  //VecAXPY(dg_new, dgO_norm, BHD->uc[1]);
	  VecAXPY(dg_new, 1, BHD->uc[1]);
/* 	  dgO_norm = ImposeBoundaryConditionII(BHD, dg_new, zpad); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " %e ", dgO_norm); */
	  //VecShift(dg_new, -dgO_norm);
	  //ImposeBoundaryCondition( BHD, dg_new);
	  ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xO, zpad, NULL);
	  Zeropad_Function(BHD, dg_new, zpad, 0.0);
	  //Smooth_Function(BHD, dg_new, zpad-1, zpad, 0.0);

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);  */
/* 	  exit(1);  */

	  VecCopy(dg_new, dg_newO);


/* 	  VecCopy(dgO, f); */
/* 	  VecAXPBY(dgO, a, (1-a), dg_new); */
/* 	  VecAXPY(f, -1.0,  dgO); */
/* 	  VecNorm(f, NORM_INFINITY, &dgO_norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"O= %e  ", dgO_norm/a); */
	  //ComputeH2O_g( gO, g0O, dgO);

/* 	  if(iter<=10) */
/* 	    a=0.01; */
/* 	  else */
/* 	    a=0.05; */
/* 	  if( iter==0) */
/* 	    { */
/* 	      VecScale(dg_newO, a); */
/* 	      VecScale(dg_newH, a); */
/* 	      EnforceNormalizationCondition(BHD, dg_newO, dg_newH, gO, gH); */
/* 	      a=1.0; */
/* 	    } */

	  /* Move dgH */
	  VecCopy(dgH, f);
/* 	  VecAXPY(dg_histH, 1.0, dg_newH); */
/* 	  norm = (iter)%50 + 1; */

 	  VecAXPBY(dgH, a, (1-a), dg_newH);
 	  VecAXPY(f, -1.0, dgH);
 	  VecNorm(f, NORM_INFINITY, &dgH_norm);

 	  PetscPrintf(PETSC_COMM_WORLD,"H= %e (a=%f) ", dgH_norm/a, a);

	  /* Move dgO */
	  if(0&& iter>10)
	    {
	      aO = GetOptimalStepSize(BHD, dg_newO, dgO, dgH);
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
/* 	      VecAXPY(dg_histO, 1.0, dg_newO); */
/* 	      norm = (iter)%20 + 1; */
	      VecAXPBY(dgO, a, (1-a), dg_newO);
	      VecAXPY(f, -1.0,  dgO);
	      VecNorm(f, NORM_INFINITY, &dgO_norm);
	      PetscPrintf(PETSC_COMM_WORLD,"O= %e (a=%f) ", dgO_norm/a, a);
	    }
	  //a=a+0.002;
/* 	  if( !((iter+1)%50) && iter>0) */
/* 	    { */
/* 	      a=0.01; */
/* 	      PetscPrintf(PETSC_COMM_WORLD,"Adding history_vec..."); */
/* 	      VecAXPY(dgO, -2.0, dg_histO); */
/* 	      VecAXPY(dgH, -2.0, dg_histH); */
/* 	      VecSet(dg_histO, 0.0); */
/* 	      VecSet(dg_histH, 0.0); */
/* 	    } */
	  ComputeH2O_g( gH, g0H, dgH);
	  ComputeH2O_g( gO, g0O, dgO);
	  norm = ComputeCharge(BHD, gH, gO);
	  PetscPrintf(PETSC_COMM_WORLD, " %e ", norm);
	  //EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH);


/* 	  if( !(iter%10) &&iter>0  ) */
/* 	    { */
/* 	      EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH); */

/* 	    } */

/* 	  if(dgO_norm/a<0.1 && dgH_norm/a<0.1) */
/* 	    { */
/* 	      EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH); */
/* 	      a0=0.1; */
/* 	    } */

/* 	  if( fabs(norm-1.0) > 0.05 ) */
/* 	    a=0.005; */
/* 	  else  */
/* 	    a=0.01; */
	   /* (fancy) step size control */

          /*
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
             (dgH_old < dgH_norm / a || dgO_old < dgO_norm / a))
	    {
	      upwards = 1;
	    }
	  else if (iter > 20 && !((iter - 1) % 10) && upwards == 0 &&
		  (dgH_old < dgH_norm / a || dgO_old < dgO_norm / a))
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
	  PetscPrintf(PETSC_COMM_WORLD,"count= %d  upwards= %d", mycount, upwards);
	  dgH_old = dgH_norm / a;
	  dgO_old = dgO_norm / a;

	  /*********************************/

	  PetscPrintf(PETSC_COMM_WORLD,"\n");

	  if(dgH_norm/a<=norm_tol &&  dgO_norm/a<=norm_tol ) //&& NORM_REG<5.0e-2)
	    break;

	}
      /*************************************/
      /* output */
      namecount++;
      sprintf(nameH, "vecH-%d.m", namecount-1);
      sprintf(nameO, "vecO-%d.m", namecount-1);

      PetscPrintf(PETSC_COMM_WORLD,"Writing files...");
      bgy3d_save_vec_ascii (nameH, gH); /* g_H */
      bgy3d_save_vec_ascii (nameO, gO); /* g_O */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
      /************************************/

      /* save g to binary file */
      if (bgy3d_getopt_test ("-saveH2O")) {
	  PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
          bgy3d_save_vec ("dgH.bin", dgH); /* dgH */
          bgy3d_save_vec ("dgO.bin", dgO); /* dgO */
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

  BGY3dH2OData_free2(BHD);


  //VecView(BHD->uc[1],PETSC_VIEWER_STDERR_WORLD);
   //exit(1);

  return PETSC_NULL;
}



Vec BGY3dM_solve_H2O_3site(ProblemData *PD, Vec g_ini, int vdim)
{
  State *BHD;
  real a0=0.1, a1, a, damp_start=0.0, norm_tol=1.0e-2, zpad=1000.0, damp, damp_LJ;
  real count=0.0, norm, aO;
  int max_iter=100, iter;
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

  BHD = BGY3dH2OData_malloc(PD);
  BHD->rho_H = 2.*BHD->rho_H;
  if( r_HH <0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Solvent not a 3-Site model!\n");
      exit(1);
    }

  /* read BGY3d specific things from command line */
  /* Mixing parameter */
  bgy3d_getopt_real ("-lambda", &a0);
  if(a0>1 || a0<0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lambda out of range: lambda=%f\n",a0);
      exit(1);
    }
   /* Get damp_start from command line*/
  bgy3d_getopt_real ("-damp_start", &damp_start);
  /* Number of total iterations */
  bgy3d_getopt_int ("-max_iter", &max_iter);
  /* norm_tol for convergence test */
  bgy3d_getopt_real ("-norm_tol", &norm_tol);
  /* Zeropad */
  bgy3d_getopt_real ("-zpad", &zpad);
  /*********************************/

  PetscPrintf(PETSC_COMM_WORLD,"lambda = %f\n",a0);
  PetscPrintf(PETSC_COMM_WORLD,"tolerance = %e\n",norm_tol);
  PetscPrintf(PETSC_COMM_WORLD,"zpad = %f\n",zpad);
  PetscPrintf(PETSC_COMM_WORLD,"max_iter = %d\n",max_iter);

  ImposeBoundaryCondition_Initialize( BHD, zpad);
#ifdef L_BOUNDARY
  BHD->zpad = zpad;
  /* Assemble Laplacian matrix */
  InitializeLaplaceMatrix(BHD, zpad);
  /* Create KSP environment */
  InitializeKSPSolver(BHD);
#endif
#ifdef L_BOUNDARY_MG
  BHD->zpad = zpad;
  /* Create KSP environment */
  InitializeDMMGSolver(BHD);
#endif

  DACreateGlobalVector(BHD->da, &gH);
  DACreateGlobalVector(BHD->da, &gO);
  DACreateGlobalVector(BHD->da, &dgH);
  DACreateGlobalVector(BHD->da, &dgO);
  DACreateGlobalVector(BHD->da, &dg_new);
  DACreateGlobalVector(BHD->da, &dg_new2);
  DACreateGlobalVector(BHD->da, &f);

  DACreateGlobalVector(BHD->da, &tH);
  DACreateGlobalVector(BHD->da, &tO);

  DACreateGlobalVector(BHD->da, &dg_newH);
  DACreateGlobalVector(BHD->da, &dg_newO);

  DACreateGlobalVector(BHD->da, &dg_histO);
  DACreateGlobalVector(BHD->da, &dg_histH);
  VecSet(dg_histH, 0.0);
  VecSet(dg_histO, 0.0);

  g0H=BHD->g_ini[0];
  g0O=BHD->g_ini[1];

/*   VecView(BHD->fHO_l[0],PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);   */

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dg_new,0.0);

  if (bgy3d_getopt_test ("-fromg2")) {
      ComputedgFromg (dgH, g0H, BHD->g2HO);
      ComputedgFromg (dgO, g0O, BHD->g2O);
  }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("-loadH2O")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      bgy3d_load_vec ("dgH.bin", &dgH); /* dgH */
      bgy3d_load_vec ("dgO.bin", &dgO); /* dgO */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  }

  for( damp=damp_start; damp <=1; damp+=0.1)
    {
      if(damp==-0.01)
	{
	  damp_LJ=0;
	  //a0=0.4;
	  RecomputeInitialFFTs(BHD, 0.0, 1.0);
	  //RecomputeInitialSoluteData_Water(BHD, 0.0, 1.0, zpad);
	  //RecomputeInitialSoluteData_Hexane(BHD, 0.0, 1.0, zpad);
	  //RecomputeInitialSoluteData_Methanol(BHD, 0.0, 1.0, zpad);
	  RecomputeInitialSoluteData_ButanoicAcid(BHD, 0.0, 1.0, zpad);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}
      else if(damp==0.0)
	{
	  damp_LJ=1.0;
	  //a0=0.5;
	  RecomputeInitialFFTs(BHD, 0.0, 1.0);
	  //RecomputeInitialSoluteData_Water(BHD, 0.0, 1.0, zpad);
	  //RecomputeInitialSoluteData_Hexane(BHD, 0.0, 1.0, zpad);
	  //RecomputeInitialSoluteData_Methanol(BHD, 0.0, 1.0, zpad);
	  RecomputeInitialSoluteData_ButanoicAcid(BHD, 0.0, 1.0, zpad);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}
      else
	{
	  damp_LJ=1.0;
	  //a0=0.0002/damp;
	  count+=1.0;
	  //a0 = 0.1/4./(count);
	  a0 = 0.1/(count+5.0);
	  RecomputeInitialFFTs(BHD, (damp), 1.0);
	  //RecomputeInitialSoluteData_Water(BHD, (damp), 1.0, zpad);
	  //RecomputeInitialSoluteData_Hexane(BHD, (damp), 1.0, zpad);
	  //RecomputeInitialSoluteData_Methanol(BHD, (damp), 1.0, zpad);
	  RecomputeInitialSoluteData_ButanoicAcid(BHD, (damp), 1.0, zpad);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}

      //Smooth_Function(BHD, g0O, zpad-1, zpad, 0.0);
      //Smooth_Function(BHD, g0H, zpad-1, zpad, 0.0);

      ImposeLaplaceBoundary(BHD, g0H, tH, BHD->xH, zpad, NULL);
      ImposeLaplaceBoundary(BHD, g0O, tH, BHD->xO, zpad, NULL);
      Zeropad_Function(BHD, g0O, zpad, 0.0);
      Zeropad_Function(BHD, g0H, zpad, 0.0);
      /* g=g0*exp(-dg) */

      ComputeH2O_g( gH, g0H, dgH);
      ComputeH2O_g( gO, g0O, dgO);


  /*     VecCopy(BHD->g2HO, gH);  */
/*       VecCopy(BHD->g2O, gO);  */

      a=a0;
      a1=a0;
      for(iter=0; iter<max_iter; iter++)
	{

	  if( !(iter%10) && iter>0 )
	    a=a1;
	  else
	    a=a0;

/* 	  if( !(iter%50) && iter>0 && NORM_REG>=5.0e-2) */
/* 	    NORM_REG/=2.; */

	  PetscPrintf(PETSC_COMM_WORLD,"iter %d: function norms: %e ", iter+1, NORM_REG);


	  /* H */
	  VecSet(dg_new,0.0);
	  Compute_H2O_interS(BHD,
			     BHD->fg2HO_fft, gO, BHD->ucHO_fft, BHD->fH_fft,
			     1.0, BHD->rho_O, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS(BHD,
			     BHD->fg2HH_fft, gH, BHD->ucH_fft, BHD->fH_fft,
			     0.0, BHD->rho_H, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  VecScale(dg_new,damp_LJ);

	  /* Coulomb long */
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2HOl_fft, gO, BHD->ucHO_fft, BHD->fH_fft,
			       1.0, BHD->rho_O, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2HHl_fft, gH, BHD->ucH_fft, BHD->fH_fft,
			       0.0, BHD->rho_H, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);

	  Solve_NormalizationH2O_smallII( BHD, gH, r_HH, gH, tH , dg_new2, f, zpad);
	  //Compute_dg_H2O_intra_lnIII(BHD, gH, tH, r_HH, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gH, tH, r_HH, dg_new2, f);
	  Compute_dg_H2O_intra_ln(BHD, tH, r_HH, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gO, tO , dg_new2, f, zpad);
	  //Compute_dg_H2O_intra_lnIII(BHD, gO, tO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gO, tO, r_HO, dg_new2, f);
	  Compute_dg_H2O_intra_ln(BHD, tO, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);

/*  	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);      */
/*  	  exit(1);    */

	  //dgH_norm = ComputeCharge(BHD, gH, gO);
	  //PetscPrintf(PETSC_COMM_WORLD, " %e ", dgH_norm);
/* 	  if(dgH_norm<1.0) */
/* 	    VecScale(dg_new, 1.2); */
/* 	  else if(dgH_norm >1.0) */
/* 	    VecScale(dg_new, 0.8); */


	  //VecShift(dg_new, -0.01);
	  //VecAXPY(dg_new, dgH_norm,  BHD->uc[0]);
	  VecAXPY(dg_new, 1.0, BHD->uc[0]);
/* 	  dgH_norm = ImposeBoundaryConditionII(BHD, dg_new, zpad); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " %e ", dgH_norm); */
	  //VecShift(dg_new, -dgH_norm);
	  //ImposeBoundaryCondition( BHD, dg_new);

	  ti=ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xH, zpad, &iteri);
	  Zeropad_Function(BHD, dg_new, zpad, 0.0);
	  //Smooth_Function(BHD, dg_new, zpad-1, zpad, 0.0);

	  PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);



	  VecCopy(dg_new, dg_newH);

/* 	  VecCopy(dgH, f); */
/* 	  VecAXPBY(dgH, a, (1-a), dg_new); */
/* 	  VecAXPY(f, -1.0, dgH); */
/* 	  VecNorm(f, NORM_INFINITY, &dgH_norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"H= %e  ", dgH_norm/a); */
	  //ComputeH2O_g( gH, g0H, dgH);


	  /* O */
	  VecSet(dg_new,0.0);
	  Compute_H2O_interS(BHD,
			     BHD->fg2OO_fft, gO, BHD->ucO_fft, BHD->fO_fft,
			     1.0, BHD->rho_O, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS(BHD,
			     BHD->fg2HO_fft, gH, BHD->ucHO_fft, BHD->fO_fft,
			     0.0, BHD->rho_H, dg_new2);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  VecScale(dg_new,damp_LJ);

	  /* Coulomb long */
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2OOl_fft, gO, BHD->ucO_fft, BHD->fO_fft,
			       1.0, BHD->rho_O, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Compute_H2O_interS_C(BHD,
			       BHD->fg2HOl_fft, gH, BHD->ucHO_fft, BHD->fO_fft,
			       0.0, BHD->rho_H, dg_new2, damp);
	  VecAXPY(dg_new, 1.0, dg_new2);



	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gH, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gH, tH, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 2.0, dg_new2);

/*  	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);      */
/*  	  exit(1);    */


	  //dgO_norm = ComputeCharge(BHD, gH, gO);
	  //PetscPrintf(PETSC_COMM_WORLD, " %e ", dgO_norm);
/* 	  if(dgO_norm<1.0) */
/* 	    VecScale(dg_new, 0.9); */
/* 	  else if(dgO_norm >1.0) */
/* 	    VecScale(dg_new, 1.1); */


	  //VecShift(dg_new, 0.01);
	  //VecAXPY(dg_new, dgO_norm, BHD->uc[1]);
	  VecAXPY(dg_new, 1, BHD->uc[1]);
/* 	  dgO_norm = ImposeBoundaryConditionII(BHD, dg_new, zpad); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " %e ", dgO_norm); */
	  //VecShift(dg_new, -dgO_norm);
	  //ImposeBoundaryCondition( BHD, dg_new);

	  //ti=ImposeLaplaceBoundary(BHD, dg_new, tH, tO, zpad, &iteri);

	  ti=ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xO, zpad, &iteri);
	  Zeropad_Function(BHD, dg_new, zpad, 0.0);
	  //Smooth_Function(BHD, dg_new, zpad-1, zpad, 0.0);
	  PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);

/* 	  if(iter==51) */
/* 	    { */
/* 	      VecView(tO,PETSC_VIEWER_STDERR_WORLD);  */
/* 	      exit(1); */
/* 	    }  */

	  VecCopy(dg_new, dg_newO);


/* 	  VecCopy(dgO, f); */
/* 	  VecAXPBY(dgO, a, (1-a), dg_new); */
/* 	  VecAXPY(f, -1.0,  dgO); */
/* 	  VecNorm(f, NORM_INFINITY, &dgO_norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"O= %e  ", dgO_norm/a); */
	  //ComputeH2O_g( gO, g0O, dgO);

/* 	  if(iter<=10) */
/* 	    a=0.01; */
/* 	  else */
/* 	    a=0.05; */
/* 	  if( iter==0) */
/* 	    { */
/* 	      VecScale(dg_newO, a); */
/* 	      VecScale(dg_newH, a); */
/* 	      EnforceNormalizationCondition(BHD, dg_newO, dg_newH, gO, gH); */
/* 	      a=1.0; */
/* 	    } */

	  /* Move dgH */
	  VecCopy(dgH, f);
/* 	  VecAXPY(dg_histH, 1.0, dg_newH); */
/* 	  norm = (iter)%50 + 1; */

 	  VecAXPBY(dgH, a, (1-a), dg_newH);
 	  VecAXPY(f, -1.0, dgH);
 	  VecNorm(f, NORM_INFINITY, &dgH_norm);

 	  PetscPrintf(PETSC_COMM_WORLD,"H= %e (a=%f) ", dgH_norm/a, a);

	  /* Move dgO */
	  if(0&& iter>10)
	    {
	      aO = GetOptimalStepSize(BHD, dg_newO, dgO, dgH);
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
/* 	      VecAXPY(dg_histO, 1.0, dg_newO); */
/* 	      norm = (iter)%20 + 1; */
	      VecAXPBY(dgO, a, (1-a), dg_newO);
	      VecAXPY(f, -1.0,  dgO);
	      VecNorm(f, NORM_INFINITY, &dgO_norm);
	      PetscPrintf(PETSC_COMM_WORLD,"O= %e (a=%f) ", dgO_norm/a, a);
	    }
	  //a=a+0.002;
/* 	  if( !((iter+1)%50) && iter>0) */
/* 	    { */
/* 	      a=0.01; */
/* 	      PetscPrintf(PETSC_COMM_WORLD,"Adding history_vec..."); */
/* 	      VecAXPY(dgO, -2.0, dg_histO); */
/* 	      VecAXPY(dgH, -2.0, dg_histH); */
/* 	      VecSet(dg_histO, 0.0); */
/* 	      VecSet(dg_histH, 0.0); */
/* 	    } */
	  ComputeH2O_g( gH, g0H, dgH);
	  ComputeH2O_g( gO, g0O, dgO);
	  norm = ComputeCharge(BHD, gH, gO);
	  PetscPrintf(PETSC_COMM_WORLD, " %e ", norm);
	  //EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH);


/* 	  if( !(iter%10) &&iter>0  ) */
/* 	    { */
/* 	      EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH); */

/* 	    } */

/* 	  if(dgO_norm/a<0.1 && dgH_norm/a<0.1) */
/* 	    { */
/* 	      EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH); */
/* 	      a0=0.1; */
/* 	    } */

/* 	  if( fabs(norm-1.0) > 0.05 ) */
/* 	    a=0.005; */
/* 	  else  */
/* 	    a=0.01; */
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
      sprintf(nameH, "vecH-%d.m", namecount-1);
      sprintf(nameO, "vecO-%d.m", namecount-1);

      PetscPrintf(PETSC_COMM_WORLD,"Writing files...");
      bgy3d_save_vec_ascii (nameH, gH); /* g_H */
      bgy3d_save_vec_ascii (nameO, gO); /* g_O */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");

      /************************************/

      /* save g to binary file */
      if (bgy3d_getopt_test ("-saveH2O")) {
	  PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
          bgy3d_save_vec ("dgH.bin", dgH); /* dgH */
          bgy3d_save_vec ("dgO.bin", dgH); /* dgO */
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

  BGY3dH2OData_free2(BHD);


  //VecView(BHD->uc[1],PETSC_VIEWER_STDERR_WORLD);
   //exit(1);

  return PETSC_NULL;
}
