/*==========================================================*/
/*  $Id: bgy3d.c,v 1.48 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/


/*==========================================*/
/* Hinweise zur Faltung: */
/* g(x) muss immer um 0 herum bekannt sein. Also sollten auch */
/* die x_M um 0 herum liegen. */
/*===============================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3dH2OS.h"          /* BGY3d_solve_H2O_2site */
#include "bgy3dH2O.h"           /* BGY3d_solve_2site */
#include "bgy3dmolecule.h"      /* BGY3d_solve_DiatomicAB */
#include "bgy3dH2ONewton.h"

static char helptext[] = "Solving BGY3d equation.\n";

int verbosity=0;

typedef Vec Solver(ProblemData *PD, Vec g_ini, int vdim);

int main(int argc, char **argv)
{
  ProblemData PD;
  int ierr, N=0;
  real beta=0.6061, rho=0.3, h=0.5, interval[2]={-5.0,5.0};
  real mpi_start, mpi_stop;
  int np;
  Vec g, g_ini;
  Solver *s=NULL;

  verbosity=0;

  PetscInitialize( &argc, &argv, (char*)0, helptext);


  MPI_Comm_size(PETSC_COMM_WORLD, &np);
  PetscPrintf(PETSC_COMM_WORLD, "NP= %d\n", np);


  /* computation time measurement start point */
  MPI_Barrier( PETSC_COMM_WORLD);
  mpi_start = MPI_Wtime();

  /* Read the command  line options. Petsc insists on  keys having the
     leading dash, so keep them for the moment. */

  /* Grid points in 1 dimension */
  bgy3d_getopt_int ("-N", &N);

  /* inverse temperature */
  bgy3d_getopt_real ("-beta", &beta);

  /* Density */
  bgy3d_getopt_real ("-rho", &rho);

  /* set verbosity */
  bgy3d_getopt_int ("-v", &verbosity);

  /* mesh width */
  bgy3d_getopt_real ("-mesh", &h);

  /*=================================*/
  /* set Problem Data */
  /*=================================*/

  if(N==0)
    N= (int) ceil((interval[1]-interval[0])/h);
  else
    h = (interval[1]-interval[0])/N;
  FOR_DIM
    PD.N[dim] = N;
  PD.N3 = PD.N[0]*PD.N[1]*PD.N[2];
  FOR_DIM
    PD.h[dim] = h;
  PD.interval[0] = interval[0];
  PD.interval[1] = interval[1];
  PD.beta = beta;
  PD.rho  = rho;
  PD.g_xm = 1.0; /* deprecated */

  PData_CreateParallel(&PD);
  PetscViewerSetFormat(PETSC_VIEWER_STDERR_WORLD,PETSC_VIEWER_ASCII_MATLAB);
  PetscViewerSetFormat(PETSC_VIEWER_STDERR_SELF,PETSC_VIEWER_ASCII_MATLAB);
  /*==================================*/

  PetscPrintf(PETSC_COMM_WORLD, "Grid size N=%d %d %d\n",PD.N[0], PD.N[1], PD.N[2]);
  PetscPrintf(PETSC_COMM_WORLD, "Total dof N^3=%d\n",PD.N3);
  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", interval[0], interval[1]);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", h);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", beta);
  PetscPrintf(PETSC_COMM_WORLD, "rho = %f\n", rho);

  //PetscPrintf(PETSC_COMM_WORLD, "\tATTENTION: Factor 2 is included!!! But why???\n");

  //if(PD.id==1)
/*   start_debugger(); */
/*   sleep(5); */


  /* Read method to solve from command line */
  if (bgy3d_getopt_test ("-simple"))
    s = BGY3d_solve;

  /* if (bgy3d_getopt_test ("-full")) */
  /*     s = BGY3d_vec_solve; */

  if (bgy3d_getopt_test ("-HNC"))
    s = HNC3d_Solve_h;

  if (bgy3d_getopt_test ("-HNCNewton"))
    s = HNC3dNewton2_solve;

  if (bgy3d_getopt_test ("-DIV"))
    s = BGY3dDiv_solve2;

  if (bgy3d_getopt_test ("-DIV2"))
    s = BGY3dDiv_solve;

  if (bgy3d_getopt_test ("-BGYTEST"))
    s = BGY3dDiv_test;

  if (bgy3d_getopt_test ("-BGYFOURIER"))
    s =  BGY3dDiv_solve_Fourier;

  if (bgy3d_getopt_test ("-BGYFOURIERTEST"))
    s =  BGY3dDiv_solve_FourierTest;

  if (bgy3d_getopt_test ("-BGYCONVOLUTIONTEST"))
    s =  BGY3d_Convolution_Test;

  if (bgy3d_getopt_test ("-BGYDIATOMIC"))
    s =  BGY3d_solve_DiatomicAB;

  if (bgy3d_getopt_test ("-BGY2Site"))
    s =  BGY3d_solve_2site;

  if (bgy3d_getopt_test ("-BGY3Site"))
    s =  BGY3d_solve_3site;

  if (bgy3d_getopt_test ("-BGYM2Site"))
    s =  BGY3dM_solve_H2O_2site;

  if (bgy3d_getopt_test ("-BGYM3Site"))
    s =  BGY3dM_solve_H2O_3site;

  if (bgy3d_getopt_test ("-BGYH2ONEWTON"))
    s = BGY3d_SolveNewton_H2O;

  if (bgy3d_getopt_test ("-BGYH2OSNEWTON"))
    s = BGY3d_SolveNewton_H2OS;

  if (bgy3d_getopt_test ("-BGYH2OSFNEWTON"))
    s = BGY3d_SolveNewton_H2OSF;

  if(s)
    {
      /* load initial configuration from file ??? */
      if (bgy3d_getopt_test ("-load")) {
          bgy3d_load_vec ("g.bin", &g_ini);
          PetscPrintf(PETSC_COMM_WORLD,"g_ini loaded from file \"g.bin\".\n");

          g= (*s)(&PD, g_ini, 0);
      }
      else
        g= (*s)(&PD, PETSC_NULL, 0);

      /* computation time measurement end point*/
      MPI_Barrier( PETSC_COMM_WORLD);
      mpi_stop = MPI_Wtime();
      PetscPrintf(PETSC_COMM_WORLD,"Total computation time: %.4f s\n",
                  mpi_stop-mpi_start);



      /* Output result */
      if( g != PETSC_NULL) {
          bgy3d_save_vec_ascii ("vec.m", g);

          /* save g to binary file */
          if (bgy3d_getopt_test ("-save")) {
              bgy3d_save_vec ("g.bin", g);
	      PetscPrintf(PETSC_COMM_WORLD,"Result written to file \"g.bin\".\n");
          }

	  VecDestroy(g);
      }

    }
  else
    PetscPrintf(PETSC_COMM_WORLD, "Please choose one of: -BGYFOURIER, -BGYTEST, -BGYDIATOMIC, -HNC or -HNCNewton !\n");

  ierr = PetscFinalize();CHKERRQ(ierr);

  /*
   * Original version of BGY3D  executable returned 1. Make interprets
   * non-zero codes as a failure, so this was changed to 0:
   */
  return 0;
}



// real Lennard_Jones(real r, void *LJ_params)
real Lennard_Jones(real r, real epsilon, real sigma)
{
  // real epsilon, sigma, sr6, sr, re;
  real sr6, sr, re;

  r=r+SHIFT;

  // epsilon = ((double*)LJ_params)[0];
  // sigma   = ((double*)LJ_params)[1];

  sr = sigma/r;
  sr6 = SQR(sr)*SQR(sr)*SQR(sr);

  re= 4.*epsilon*sr6*(sr6-1.);

  if(fabs(re)>epsilon*CUTOFF)
    return epsilon*CUTOFF;
  else
    return re;
}

// real Lennard_Jones_grad(real r, real xr, void *LJ_params)
real Lennard_Jones_grad(real r, real xr, real epsilon, real sigma)
{
  // real epsilon, sigma, sr6, sr, re;
  real sr6, sr, re;

  r=r+SHIFT;

  // epsilon = ((double*)LJ_params)[0];
  // sigma   = ((double*)LJ_params)[1];

  if(xr==0)
    return 0;
  if(r==0)
    return -epsilon*CUTOFF;

  sr = sigma/r;
  sr6 = SQR(sr)*SQR(sr)*SQR(sr);

  re = -24.*epsilon*sr6/r*(2.*sr6-1.)*xr/r;

  //re = -24.*epsilon*pow(sigma,6.0)/pow(r,7.0)*(2.*pow(sigma,6.0)/pow(r,6.0)-1.)*xr/r;

  if(re>epsilon*CUTOFF)
    return epsilon*CUTOFF;
  else if(re<-epsilon*CUTOFF)
    return -epsilon*CUTOFF;
  else
    return re;
}

void PData_CreateParallel(ProblemData *PD)
{

  MPI_Comm_size(PETSC_COMM_WORLD, &(PD->np));
  MPI_Comm_rank(PETSC_COMM_WORLD, &(PD->id));
}


BGY3dParameter BGY3dParameter_malloc(ProblemData *PD, int vdim)
{
  BGY3dParameter params;
  DA da;
  int n[3], x[3], i[3], N_g2, index, N_M, k, ic1, ic2, ic3, N[3];
  PetscScalar ***x_vec, ***force_vec, ***Ftimesg2_vec, ***boundary_vec,h[3],
    ***force_single_vec, interval[2];
  PetscScalar r[3], r_s, L, *g2_vec, h_g2;
  PetscViewer pview;
  Vec g2;
  real **x_M;
  int bufsize;
  // Local LJ params
  real epsilon, sigma;

  params = (BGY3dParameter) malloc(sizeof(*params));

  params->vec_dim = vdim;
  // params->LJ_params = (void* ) malloc(sizeof(real)*2);
  // ((real*)(params->LJ_params))[0] = 1.0;   /* espilon */
  // ((real*)(params->LJ_params))[1] = 1.0;   /* sigma   */
  params->LJ_params[0] = 1.0;
  params->LJ_params[1] = 1.0;
  epsilon = params->LJ_params[0];
  sigma = params->LJ_params[1];


  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];
  FOR_DIM
    h[dim]=PD->h[dim];
  FOR_DIM
    N[dim]=PD->N[dim];
  L=interval[1]-interval[0];

  /* Create distributed array */
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
	     1,1, PETSC_NULL,PETSC_NULL,PETSC_NULL, &(params->da));

  da = params->da;

  /* Create global vectors */
  DACreateGlobalVector(da, &(params->x));
  VecDuplicate(params->x, &(params->force));
  VecDuplicate(params->x, &(params->force_single));
  VecDuplicate(params->x, &(params->Ftimesg2));
  VecDuplicate(params->x, &(params->boundary));
  VecDuplicate(params->x, &(params->v1));
  VecDuplicate(params->x, &(params->v2));
  VecDuplicate(params->x, &(params->v3));
  VecDuplicate(params->x, &(params->pre));

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if( verbosity >2)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Subgrids on processes:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "id %d of %d: %d %d %d\t%d %d %d\n",
			      PD->id, PD->np, x[0], x[1], x[2], n[0], n[1], n[2]);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }

  /* Create Matrix with appropriate non-zero structure */
  DAGetMatrix( da, MATMPIAIJ, &(params->M));
  ComputeMatrixStencil(PD, da, params->M, vdim);

  //MatView(params->M,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);



  /* Load g2 from file, sequential only */
  /* g2 has to be on a grid [0,L] with L=interval[1]-interval[0] */
  PetscViewerBinaryOpen(PETSC_COMM_SELF,"g2file.bin" , FILE_MODE_READ,
			&pview);
  VecLoad( pview, VECSEQ, &g2);
  VecGetSize(g2,&N_g2);
  h_g2 = L/N_g2;


  /* Load molecule from file */
  x_M = Load_Molecule(&N_M);



  VecSet(params->force,0.0);
  DAVecGetArray(da, params->x, &x_vec);
  DAVecGetArray(da, params->force, &force_vec);
  DAVecGetArray(da, params->force_single, &force_single_vec);
  DAVecGetArray(da, params->Ftimesg2, &Ftimesg2_vec);
  DAVecGetArray(da, params->boundary, &boundary_vec);
  VecGetArray(g2, &g2_vec);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set x vector */
	  x_vec[i[2]][i[1]][i[0]] = i[0]*h[0]+interval[0];



	  /* set force vector */
	  /* loop over particles and grid */
	  for(k=0; k<N_M; k++)
	    {
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim]+interval[0]-x_M[k][dim];
/* 		  if( i[dim]>=N[dim]/2) */
/* 		    r[dim] -= L; */
		}

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );
	      force_vec[i[2]][i[1]][i[0]] +=
		// Lennard_Jones_grad( r_s, r[vdim], params->LJ_params);
		Lennard_Jones_grad( r_s, r[vdim], epsilon, sigma);
	    }

	  /* set force*g2 and force_single vector */
	  /* these vectors live on a different grid -> [0,L]^3 for fft */
	  FOR_DIM
	    {
	      r[dim] = i[dim]*h[dim];
	      if( i[dim]>=N[dim]/2)
		r[dim] -= L;
	    }

	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  force_single_vec[i[2]][i[1]][i[0]] =
	    // Lennard_Jones_grad( r_s, r[vdim], params->LJ_params);
	    Lennard_Jones_grad( r_s, r[vdim], epsilon, sigma);

	  index =(int) floor(r_s/h_g2);
	  /* g2=1 if r_s>l */
	  if(index >= N_g2-1)
	    Ftimesg2_vec[i[2]][i[1]][i[0]] =
	      // Lennard_Jones_grad( r_s, r[vdim], params->LJ_params);
	      Lennard_Jones_grad( r_s, r[vdim], epsilon, sigma);
	  else
	    {
	      Ftimesg2_vec[i[2]][i[1]][i[0]] =
		// Lennard_Jones_grad( r_s, r[vdim], params->LJ_params)*
		Lennard_Jones_grad( r_s, r[vdim], epsilon, sigma)*
		(g2_vec[index]+ (g2_vec[index+1]-g2_vec[index])/h_g2*
		 (r_s-index*h_g2));

	    }

	}



  /* set boundary vector: r:g=1 l:g=-1 */
#ifdef CENTRALDIFF
  ic1=(vdim)%3;
  ic2=(vdim+1)%3;
  ic3=(vdim+2)%3;
  if( x[ic1] == 0 )
    {
      i[ic1] = 0;
      for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	  boundary_vec[i[2]][i[1]][i[0]] = -0.5;
    }
  if( x[ic1]+n[ic1] == PD->N[ic1])
    {
      i[ic1] = PD->N[ic1]-1;
      for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	  boundary_vec[i[2]][i[1]][i[0]] = 0.5;
    }
#else
  ic1=(vdim)%3;
  ic2=(vdim+1)%3;
  ic3=(vdim+2)%3;
  if( x[ic1]+n[ic1] == PD->N[ic1])
    {
      i[ic1] = PD->N[ic1]-1;
      for(i[ic2]= x[ic2]; i[ic2]<x[ic2]+n[ic2]; i[ic2]++)
	for(i[ic3]=x[ic3] ; i[ic3]<x[ic3]+n[ic3]; i[ic3]++)
	  boundary_vec[i[2]][i[1]][i[0]] = 1;
    }
#endif




  DAVecRestoreArray(da, params->x, &x_vec);
  DAVecRestoreArray(da, params->force, &force_vec);
  DAVecRestoreArray(da, params->force_single, &force_single_vec);
  DAVecRestoreArray(da, params->Ftimesg2, &Ftimesg2_vec);
  DAVecRestoreArray(da, params->boundary, &boundary_vec);
  VecRestoreArray(g2, &g2_vec);

/*   VecView(params->boundary,PETSC_VIEWER_STDERR_WORLD); */
/*   VecView(params->force,PETSC_VIEWER_STDERR_WORLD);  */
/*   VecView(params->Ftimesg2,PETSC_VIEWER_STDERR_WORLD); */
/*    exit(1);   */



  /* Create plan for 3d fft */
  params->fft_plan = fft_3d_create_plan(PETSC_COMM_WORLD,
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
  if(params->fft_plan == NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Failed to get fft_plan of proc %d.\n",
		  PD->id);
      exit(1);
    }


  /* Compute fft for f*g2 */
  params->Ftimesg2_fft = NULL;
  params->Ftimesg2_fft = ComputeFFTfromVec(da, params->fft_plan,
  					   params->Ftimesg2,
  					   params->Ftimesg2_fft);

  params->PD=PD;

  /* Initialize matrix preconditioner */
#ifdef MATPRECOND
  params->MP = MatPrecond_malloc(params);
#endif
/*   MatView(params->MP->P,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

  VecDestroy(g2);

  Molecule_free(x_M, N_M);

  return params;
}

/* Destroy BGY3dParameter struct */
void BGY3dParameter_free(BGY3dParameter params)
{
  DADestroy(params->da);
  MatDestroy(params->M);
  VecDestroy(params->x);
  VecDestroy(params->force);
  VecDestroy(params->force_single);
  VecDestroy(params->Ftimesg2);
  VecDestroy(params->boundary);
  VecDestroy(params->v1);
  VecDestroy(params->v2);
  VecDestroy(params->v3);
  VecDestroy(params->pre);

  fft_3d_destroy_plan(params->fft_plan);

#ifdef MATPRECOND
  MatPrecond_free(params->MP);
#endif

  // free(params->LJ_params);

  free(params);
}

#ifdef MATPRECOND
MatPrecond MatPrecond_malloc(BGY3dParameter params)
{
  MatPrecond MP;
  ProblemData *PD;
  PetscScalar ***force_vec;
  int i[3], x[3], n[3], vdim;
  MatStencil col,row;
  PetscScalar v;
  real beta, h[3];
  Mat P;

  MP = (MatPrecond) malloc(sizeof(*MP));

  PD=params->PD;

  DAGetCorners(params->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));



  /* copy M to P */
  DAGetMatrix( params->da, MATMPIAIJ, &(MP->P));
  MatCopy(params->M, MP->P, SAME_NONZERO_PATTERN);

  /* insert diagonal elements */
  P = MP->P;
  beta = PD->beta;
  vdim = params->vec_dim;
  FOR_DIM
    h[dim]=PD->h[dim];

  DAVecGetArray(params->da, params->force, &force_vec);
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  col.i=i[0];
	  col.j=i[1];
	  col.k=i[2];
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  v= beta*h[vdim]*force_vec[i[2]][i[1]][i[0]];

	  MatSetValuesStencil(P,1,&row,1,&col,&v,ADD_VALUES);
	}
  DAVecRestoreArray(params->da, params->force, &force_vec);
  MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);

  /* initialize linear solver for preconditioning */
  KSPCreate(PETSC_COMM_WORLD, &(MP->ksp));
  KSPGetPC(MP->ksp, &(MP->pc));
  PCSetType(MP->pc, PCBJACOBI);

  /* successive calls to KSPSolve will not recompute
     preconditioner factorization*/
  KSPSetOperators(MP->ksp, MP->P, MP->P, SAME_PRECONDITIONER);

  KSPSetFromOptions(MP->ksp);


  return MP;
}

void MatPrecond_free( MatPrecond MP)
{
  MatDestroy(MP->P);
  KSPDestroy(MP->ksp);
}

#endif

real** Load_Molecule(int *N)
{
  FILE *fp;
  real **x_M;
  int i;


  fp = fopen("molecule","r");
  if(fp==NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD,"File not found.\n");
      exit(1);
    }

  fscanf(fp,"%d",N);
  if(verbosity>2)
    PetscPrintf(PETSC_COMM_WORLD,"Reading molecule data:\n%d particles\n", *N);

  x_M= (real**) malloc(sizeof(*x_M)*(*N));


  for(i=0; i<*N; i++)
    {
      x_M[i] = (real*) malloc(sizeof(real)*3);
      if( 3!=fscanf(fp,"%lf %lf %lf\n", &x_M[i][0], &x_M[i][1], &x_M[i][2]) )
	{
	  PetscPrintf(PETSC_COMM_WORLD,"Error reading file (molecule).\n");
	  exit(1);
	}
      if(verbosity>2)
	PetscPrintf(PETSC_COMM_WORLD,"%f %f %f\n", x_M[i][0],
		    x_M[i][1], x_M[i][2]);
    }
  fclose(fp);

  return x_M;

}


void Molecule_free( real **x_M, int N_M)
{
  int index;

  for(index=0; index<N_M; index++)
    free(x_M[index]);
  free(x_M);
}


/* Initialize M-Matrix with appropriate stencil */
void ComputeMatrixStencil(ProblemData *PD, DA da, Mat M, int vdim)
{
  int x[3], n[3], i[3], N[3];
  MatStencil col[2],row;
  PetscScalar v[2];



  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  FOR_DIM
    N[dim] = PD->N[dim];

#ifdef CENTRALDIFF
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  col[0].i=i[0];
	  col[0].j=i[1];
	  col[0].k=i[2];
	  col[1].i=i[0];
	  col[1].j=i[1];
	  col[1].k=i[2];
	  switch(vdim)
	    {
	    case 0: col[0].i -=1; col[1].i +=1;break;
	    case 1: col[0].j -=1; col[1].j +=1;break;
	    case 2: col[0].k -=1; col[1].k +=1;break;
	    default: PetscPrintf(PETSC_COMM_WORLD,"VEC_DIM value out of bounds.\n"); exit(1);
	    }

	  /* values to enter */
	  v[0]=-0.5;
	  v[1]= 0.5;
	  /* Boundary */
	  if( i[vdim]== 0 )
	    MatSetValuesStencil(M,1,&row,1,col+1,v+1,INSERT_VALUES);
	  else if( i[vdim]==N[vdim]-1)
	    MatSetValuesStencil(M,1,&row,1,col,v,INSERT_VALUES);
	  else
	    MatSetValuesStencil(M,1,&row,2,col,v,INSERT_VALUES);
	}
#else
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* position in matrix */
	  row.i=i[0];
	  row.j=i[1];
	  row.k=i[2];
	  col[0].i=i[0];
	  col[0].j=i[1];
	  col[0].k=i[2];
	  col[1].i=i[0];
	  col[1].j=i[1];
	  col[1].k=i[2];
	  switch(vdim)
	    {
	    case 0: col[1].i +=1;break;
	    case 1: col[1].j +=1;break;
	    case 2: col[1].k +=1;break;
	    default: PetscPrintf(PETSC_COMM_WORLD,"VEC_DIM value out of bounds.\n"); exit(1);
	    }

	  /* values to enter */
	  v[0]=-1;
	  v[1]= 1;
	  /* Boundary */
	  if( i[vdim]==N[vdim]-1)
	    MatSetValuesStencil(M,1,&row,1,col,v,INSERT_VALUES);
	  else
	    MatSetValuesStencil(M,1,&row,2,col,v,INSERT_VALUES);
	}
#endif

  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);

}

Vec BGY3d_solve(ProblemData *PD, Vec g_ini, int vec_dim)
{
  SNES snes;
  KSP ksp;
  PC  pc;
  BGY3dParameter params;
  PetscErrorCode ierr;
  Vec g, F;
  Mat J;

  params = BGY3dParameter_malloc(PD, vec_dim);

  ConvolutionTest(params);
  exit(1);


  /* Create the snes environment */
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp, &pc);

  /* set rtol, atol, dtol, maxits */
  // KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType(snes, SNESLS);

  /*
   * CHKERRQ()  is  a  macro that  expands  to  an  if with  a  return
   * statement  in  the if-block.   It  is  not  suitable for  use  in
   * functions  returning  anything   that  but  PetscErrorCode.  This
   * function returns a Vec.
   */

  if (bgy3d_getopt_test ("-user_precond")) { /* user-defined precond */
    /* Set user defined preconditioner */
    ierr = PCSetType(pc,PCSHELL); // CHKERRQ(ierr);
    assert (!ierr);
    ierr = PCShellSetApply(pc,Compute_Preconditioner); // CHKERRQ(ierr);
    assert (!ierr);
    ierr = PCShellSetContext(pc,params); // CHKERRQ(ierr);
    assert (!ierr);
  } else
    /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
    PCSetType( pc, PCJACOBI);

#ifdef MATPRECOND
  if (bgy3d_getopt_test ("-mat_precond")) { /* user-defined precond */
    /* Set user defined preconditioner */
    ierr = PCSetType(pc,PCSHELL); // CHKERRQ(ierr);
    assert (!ierr);
    ierr = PCShellSetApply(pc,Compute_Preconditioner_Mat); // CHKERRQ(ierr);
    assert (!ierr);
    ierr = PCShellSetContext(pc,params); // CHKERRQ(ierr);
    assert (!ierr);
  }
#endif



  VecDuplicate(params->x, &F);


  //MatDuplicate(params->M, MAT_DO_NOT_COPY_VALUES, &J);
  DAGetMatrix( params->da, MATMPIAIJ, &J);

  if(g_ini == PETSC_NULL)
    {
      VecDuplicate(params->x, &g);
      CreateInitialGuessFromg2(params, g);
      //CreateInitialGuess(params, g);

    }
  else
    g=g_ini;

/*   VecSet(g, 0.0);   */
/*   Compute_F(snes, g, F, (void*) params); */
  //VecView(g,PETSC_VIEWER_STDERR_WORLD);
  VecView(g,PETSC_VIEWER_STDERR_WORLD);
  exit(1);


  //TestPreconditioner(params->MP, g,F);


  /* set Function */
  if (bgy3d_getopt_test ("-kirkwood"))
    {
      SNESSetFunction(snes, F, Compute_F_Kirkwood, params);
    }
  else
    {
      SNESSetFunction(snes, F, Compute_F, params);
      /* set Jacobian */
      SNESSetJacobian(snes, J, J, Compute_J, params);
    }

  /* Matrix free environment */
  //MatCreateSNESMF(snes, g, &J);
  //SNESSetJacobian(snes, J, J, PETSC_NULL, PETSC_NULL);

  /* runtime options will override default parameters */
  SNESSetFromOptions(snes);

  /* solve problem */
  ierr=SNESSolve(snes, PETSC_NULL, g);





  /* write out solution */
  SNESGetSolution(snes, &g);
  //VecView(g,PETSC_VIEWER_STDERR_WORLD);
  //Compute_F_Kirkwood(snes, g, F, (void*) params);
  //VecView(F,PETSC_VIEWER_STDERR_WORLD);

  VecDestroy(F);
  MatDestroy(J);
  SNESDestroy(snes);
  BGY3dParameter_free(params);


  return g;
}


void CreateInitialGuess(BGY3dParameter params, Vec g)
{
  DA da;
  ProblemData *PD;
  int x[3], i[3], n[3], N_M, k;
  PetscScalar ***g_vec, r[3], r_s;
  real **x_M, h[3], beta, interval[2];
  real epsilon, sigma;

  epsilon = params->LJ_params[0];
  sigma = params->LJ_params[1];

  da =  params->da;
  PD = params->PD;

  FOR_DIM
    h[dim] = PD->h[dim];

  beta = PD->beta;
  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* Load molecule from file */
  x_M = Load_Molecule(&N_M);

  VecSet(g, 1.0);
  DAVecGetArray(da, g, &g_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  /* loop over particles and grid */
	  for(k=0; k<N_M; k++)
	    {
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim]+interval[0]-x_M[k][dim];
/* 		  if( i[dim]>=N[dim]/2) */
/* 		    r[dim] -= L; */
		}

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	      g_vec[i[2]][i[1]][i[0]] *=
		// exp(-beta* Lennard_Jones( r_s, params->LJ_params));
		exp(-beta* Lennard_Jones( r_s, epsilon, sigma));
	    }
	}
  DAVecRestoreArray(da, g, &g_vec);

/*   VecView(g,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */
  Molecule_free(x_M, N_M);

}


void CreateInitialGuessFromg2(BGY3dParameter params, Vec g)
{
  DA da;
  ProblemData *PD;
  int x[3], i[3], n[3], N_M, k, N_g2, index;
  PetscScalar ***g_vec, *g2_vec, r[3], r_s;
  real **x_M, h[3], h_g2, L, interval[2];
  PetscViewer pview;
  Vec g2;


  da =  params->da;
  PD = params->PD;

  FOR_DIM
    h[dim] = PD->h[dim];

  L = PD->interval[1]-PD->interval[0];
  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /* Load g2 from file, sequential only */
  /* g2 has to be on a grid [0,L] with L=interval[1]-interval[0] */
  PetscViewerBinaryOpen(PETSC_COMM_SELF,"g2file.bin" , FILE_MODE_READ,
			&pview);
  VecLoad( pview, VECSEQ, &g2);
  VecGetSize(g2,&N_g2);
  h_g2 = L/N_g2;


  /* Load molecule from file */
  x_M = Load_Molecule(&N_M);

  VecSet(g, 1.0);
  DAVecGetArray(da, g, &g_vec);
  VecGetArray(g2, &g2_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  /* loop over particles and grid */
	  for(k=0; k<N_M; k++)
	    {
	      FOR_DIM
		{
		  r[dim] = i[dim]*h[dim]+interval[0]-x_M[k][dim];
/* 		  if( i[dim]>=N[dim]/2) */
/* 		    r[dim] -= L; */
		}

	      r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	      index =(int) floor(r_s/h_g2);
	      if(index >= N_g2-1)
		g_vec[i[2]][i[1]][i[0]] *=1.0;
	      else
		g_vec[i[2]][i[1]][i[0]] *=
		g2_vec[index]+ (g2_vec[index+1]-g2_vec[index])/h_g2*
		  (r_s-index*h_g2);

	    }
	}
  DAVecRestoreArray(da, g, &g_vec);
  VecRestoreArray(g2, &g2_vec);

  VecDestroy(g2);
  Molecule_free(x_M, N_M);

}



PetscErrorCode Compute_F(SNES snes, Vec g, Vec f, void *pa)
{
  ProblemData *PD;
  DA da;
  BGY3dParameter params;
  FFT_DATA *fft_data, *fft_gFg2, *Ftimesg2_fft;
  int x[3], n[3], i, max_i, vdim;
  Vec pre;


  params = (BGY3dParameter) pa;
  pre = params->pre;

  PD = params->PD;
  da = params->da;
  vdim = params->vec_dim;

  /* get size of local grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  max_i = n[0]*n[1]*n[2];

  fft_gFg2 = (FFT_DATA*) calloc(max_i,sizeof(*fft_gFg2));
  Ftimesg2_fft = params->Ftimesg2_fft;

  /* g = fft(g) */
  fft_data = NULL;
  fft_data = ComputeFFTfromVec(da, params->fft_plan, g, fft_data);

  /* fft_data = fft_data * Ftimesg2_fft */
  for(i=0; i<max_i; i++)
    {
      fft_gFg2[i].re = fft_data[i].re*Ftimesg2_fft[i].re
	-fft_data[i].im*Ftimesg2_fft[i].im;

      fft_gFg2[i].im = fft_data[i].re*Ftimesg2_fft[i].im
	+fft_data[i].im*Ftimesg2_fft[i].re;
    }

  /* f=fft^-1(fft(g)*fft(F*g2)) */
  ComputeVecfromFFT(da, params->fft_plan, f, fft_gFg2);

/*   VecView(f,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1); */
  /****************************************/
  /* Set Preconditioner                   */
  /****************************************/
#ifndef MATPRECOND
  /*Preconditioner: pre=beta*rho*h/g^M(x_M)*integral */
/*     VecAXPY(pre,PD->h[0]*PD->h[1]*PD->h[2]* */
/* 	    PD->beta*PD->rho*PD->h[vdim]/PD->g_xm/PD->N[0]/PD->N[1]/PD->N[2] */
/* 	    , f); */
  VecSet(pre, 0.0);
  /* Preconditioner: pre=pre+beta*h*force */
  VecAXPY(pre,PD->beta*PD->h[vdim], params->force);
#ifndef CENTRALDIFF
  /* Preconditioner: pre=pre-1 */
  VecShift(pre,-1.0);
#endif
  /* Preconditioner: pre=pre^(-1) */
  VecReciprocal(pre);
#endif
  /******************************************/

  /* f=f*g */
  VecPointwiseMult(f, f, g);

  /* f=f*h^3/N^3*h[dim]*beta*rho/g^M(x_M) */
  VecScale(f, PD->h[0]*PD->h[1]*PD->h[2]*
	   PD->h[vdim]*PD->rho*PD->beta/PD->g_xm
	   /PD->N[0]/PD->N[1]/PD->N[2]);

/*   VecView(f,PETSC_VIEWER_STDERR_WORLD);  */


  /* v1 = M*g+B */
  MatMultAdd(params->M, g, params->boundary, params->v1);
  //MatMult( params->M, g, params->v1);

  /* v2 = beta*h*force*g */
  VecPointwiseMult(params->v2, g, params->force);
  VecScale(params->v2, PD->beta*PD->h[vdim]);

  /* f=f+v1+v2 */
  VecAXPY(f,1.0,params->v1);
  VecAXPY(f,1.0,params->v2);

/*   VecView(f,PETSC_VIEWER_STDERR_WORLD);    */
/*   VecView(params->v2,PETSC_VIEWER_STDERR_WORLD);  */
/*   VecView(params->force,PETSC_VIEWER_STDERR_WORLD); */
/*   VecView(params->x,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);   */

  free(fft_data);
  free(fft_gFg2);

  return 0;
}


PetscErrorCode Compute_J(SNES snes, Vec g, Mat *A, Mat *B, MatStructure *flag,
			 void *pa)
{
  ProblemData *PD;
  DA da;
  BGY3dParameter params;
  FFT_DATA *fft_data, *fft_gFg2, *Ftimesg2_fft;
  int x[3], n[3], i[3], max_i, vdim, k, N[3];
  Vec  f;
  PetscScalar ***f_vec, v[3];
  MatStencil col[3],row;


  row.c=0;
  col[0].c=0;
  col[1].c=0;

  *flag=SAME_NONZERO_PATTERN;

  params = (BGY3dParameter) pa;
  f = params->v3;

  PD = params->PD;
  da = params->da;
  vdim = params->vec_dim;

  FOR_DIM
    N[dim] = PD->N[dim];

  /* get size of local grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  max_i = n[0]*n[1]*n[2];

  fft_gFg2 = (FFT_DATA*) calloc(max_i,sizeof(*fft_gFg2));
  Ftimesg2_fft = params->Ftimesg2_fft;

  /* g = fft(g) */
  fft_data = NULL;
  fft_data = ComputeFFTfromVec(da, params->fft_plan, g, fft_data);

  /* fft_data = fft_data * Ftimesg2_fft */
  for(k=0; k<max_i; k++)
    {
      fft_gFg2[k].re = fft_data[k].re*Ftimesg2_fft[k].re
	-fft_data[k].im*Ftimesg2_fft[k].im;

      fft_gFg2[k].im = fft_data[k].re*Ftimesg2_fft[k].im
	+fft_data[k].im*Ftimesg2_fft[k].re;
    }

  /* f=fft^-1(fft(g)*fft(F*g2)) */
  ComputeVecfromFFT(da, params->fft_plan, f, fft_gFg2);



  /* f*=beta*rho*h/g^M(x_M)*/
  VecScale(f,PD->h[0]*PD->h[1]*PD->h[2]*
	  PD->beta*PD->rho*PD->h[vdim]/PD->g_xm/PD->N[0]/PD->N[1]/PD->N[2]);
  /* f=f+beta*h*force */
  VecAXPY(f,PD->beta*PD->h[vdim], params->force);

#ifndef CENTRALDIFF
  /*f=f-1 */
  VecShift(f,-1.0);
#endif

  /* insert entries into J */
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
	  col[0].i=i[0];
	  col[0].j=i[1];
	  col[0].k=i[2];
	  col[1].i=i[0];
	  col[1].j=i[1];
	  col[1].k=i[2];
	  col[2].i=i[0];
	  col[2].j=i[1];
	  col[2].k=i[2];
#ifdef CENTRALDIFF
	  switch(vdim)
	    {
	    case 0: col[0].i -=1;col[2].i +=1;break;
	    case 1: col[0].j -=1;col[2].j +=1;break;
	    case 2: col[0].k -=1;col[2].k +=1;break;
	    default: PetscPrintf(PETSC_COMM_WORLD,"VEC_DIM value out of bounds.\n"); exit(1);
	    }
	  /* values to enter */
	  v[0]= -0.5;
	  v[1]= f_vec[i[2]][i[1]][i[0]];
	  v[2]=  0.5;
#else
	  switch(vdim)
	    {
	    case 0: col[1].i +=1;break;
	    case 1: col[1].j +=1;break;
	    case 2: col[1].k +=1;break;
	    default: PetscPrintf(PETSC_COMM_WORLD,"VEC_DIM value out of bounds.\n"); exit(1);
	    }
	  /* values to enter */
	  v[0]= f_vec[i[2]][i[1]][i[0]];
	  v[1]= 1;
#endif



#ifdef CENTRALDIFF
	  /* Boundary */
	  if( i[vdim]== 0)
	    MatSetValuesStencil(*B,1,&row,2,col+1,v+1,INSERT_VALUES);
	  else if( i[vdim]== N[vdim]-1)
	    MatSetValuesStencil(*B,1,&row,2,col,v,INSERT_VALUES);
	  else
	    MatSetValuesStencil(*B,1,&row,3,col,v,INSERT_VALUES);
#else
	  /* Boundary */
	  if( i[vdim]==N[vdim]-1)
	    MatSetValuesStencil(*B,1,&row,1,col,v,INSERT_VALUES);
	  else
	    MatSetValuesStencil(*B,1,&row,2,col,v,INSERT_VALUES);
#endif

	}
  MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);
  DAVecRestoreArray(da, f, &f_vec);


  //VecView(params->v1,PETSC_VIEWER_STDERR_WORLD);
/*   MatView(*B,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */
  free(fft_data);
  free(fft_gFg2);

  return 0;
}


PetscErrorCode Compute_F_Kirkwood(SNES snes, Vec g, Vec f, void *pa)
{
  ProblemData *PD;
  DA da;
  BGY3dParameter params;
  FFT_DATA *fft_fg, *fft_g, *fft_fgg;
  int x[3], n[3], i, max_i, vdim;
  Vec pre, v1, v2;


  params = (BGY3dParameter) pa;
  pre = params->pre;

  PD = params->PD;
  da = params->da;
  vdim = params->vec_dim;
  v1 = params->v1;
  v2 = params->v2;


  /* get size of local grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  max_i = n[0]*n[1]*n[2];

  fft_fgg = (FFT_DATA*) calloc(max_i,sizeof(*fft_fgg));

  /* g = fft(g) */
  fft_g = NULL;
  fft_g = ComputeFFTfromVec(da, params->fft_plan, g, fft_g);
  /* f*g = fft(f*g) */
  VecPointwiseMult(v1, g, params->force_single);
  fft_fg = NULL;
  fft_fg = ComputeFFTfromVec(da, params->fft_plan, v1, fft_fg);

  /* fft(g) * fft(f*g) */
  for(i=0; i<max_i; i++)
    {
      fft_fgg[i].re = fft_g[i].re*fft_fg[i].re
	-fft_g[i].im*fft_fg[i].im;

      fft_fgg[i].im = fft_g[i].re*fft_fg[i].im
	+fft_g[i].im*fft_fg[i].re;
    }

  /* f=fft^-1(fft(g)*fft(F*g2)) */
  ComputeVecfromFFT(da, params->fft_plan, f, fft_fgg);


  /****************************************/
  /* Set Preconditioner                   */
  /****************************************/
#ifndef MATPRECOND
  /* Preconditioner: pre=beta*rho*h/g^M(x_M)*integral */
  VecAXPY(pre,PD->h[0]*PD->h[1]*PD->h[2]*
	  PD->beta*PD->rho*PD->h[vdim]/PD->N[0]/PD->N[1]/PD->N[2]
	  , f);
  /* Preconditioner: pre=pre+beta*h*force */
  VecAXPY(pre,PD->beta*PD->h[vdim], params->force);
#ifndef CENTRALDIFF
  /* Preconditioner: pre=pre-1 */
  VecShift(pre,-1.0);
#endif
  /* Preconditioner: pre=pre^(-1) */
  VecReciprocal(pre);
#endif
  /******************************************/

  /* f=f*g */
  VecPointwiseMult(f, f, g);

  /* f=f*h^3/N^3*h[dim]*beta*rho/g^M(x_M) */
  VecScale(f, PD->h[0]*PD->h[1]*PD->h[2]*
	   PD->h[vdim]*PD->rho*PD->beta
	   /PD->N[0]/PD->N[1]/PD->N[2]);

  //VecView(f,PETSC_VIEWER_STDERR_WORLD);


  /* v1 = M*g+B */
  MatMultAdd(params->M, g, params->boundary, v1);
  //MatMult( params->M, g, params->v1);

  /* v2 = beta*h*force*g */
  VecPointwiseMult(v2, g, params->force);
  VecScale(v2, PD->beta*PD->h[vdim]);

  /* f=f+v1+v2 */
  VecAXPY(f,1.0,params->v1);
  VecAXPY(f,1.0,params->v2);



/*   VecView(g,PETSC_VIEWER_STDERR_WORLD);  */
/*   VecView(f,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

  free(fft_g);
  free(fft_fg);
  free(fft_fgg);

  return 0;
}





PetscErrorCode Compute_Preconditioner(void *pa,Vec x,Vec y)
{
  BGY3dParameter params;
  PetscErrorCode ierr;

  params = (BGY3dParameter) pa;
  ierr=VecPointwiseMult(y,params->pre,x);


/*   VecView(x,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

  return ierr;
}

#ifdef MATPRECOND
PetscErrorCode Compute_Preconditioner_Mat(void *pa,Vec x,Vec y)
{
  BGY3dParameter params;
  MatPrecond MP;
  PetscErrorCode ierr;

  params = (BGY3dParameter) pa;
  MP = params->MP;

/*   VecView(x,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

  ierr = KSPSolve(MP->ksp, x, y);


  //VecView(params->pre,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);

  return ierr;
}
#endif


FFT_DATA *ComputeFFTfromVec(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft)

{
  int index=0, i[3];
  int x[3], n[3];
  PetscScalar ***g_vec;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if(g_fft==NULL)
    g_fft = (FFT_DATA*) calloc(n[0]*n[1]*n[2],sizeof(*g_fft));

  DAVecGetArray(da, g, &g_vec);
  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index=0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  //g_fft[index].re = g_vec[i[2]][i[1]][i[0]]/2.+c;
	  g_fft[index].re = g_vec[i[2]][i[1]][i[0]];
	  g_fft[index].im = 0;
	  index++;
	}
  DAVecRestoreArray(da, g, &g_vec);
  /* forward fft */
  fft_3d(g_fft, g_fft, 1, fft_plan);

  return g_fft;

}


void ComputeVecfromFFT(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft)

{
  int index=0, i[3];
  int x[3], n[3];
  PetscScalar ***g_vec;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if(g_fft==NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Error: g_fft==NULL!\n");
      exit(1);
    }

  /* backward fft */
  fft_3d(g_fft, g_fft, -1, fft_plan);

  DAVecGetArray(da, g, &g_vec);
  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index=0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  //g_vec[i[2]][i[1]][i[0]] = g_fft[index].re*2.+c; // Factor 2!!!
	  g_vec[i[2]][i[1]][i[0]] = g_fft[index].re;
	  index++;
	}
  DAVecRestoreArray(da, g, &g_vec);


}


  /* starts the gdb debugger for each process */
int start_debugger(void )
{
  pid_t pid;
  char s[100];

  pid = getpid();

  sprintf(s, "xterm -e gdb /mount/kamino/jager/work/BGY3d/test/bgy3d %d &", pid);

  return system(s);
}


void ConvolutionTest(BGY3dParameter params)
{
  ProblemData *PD;
  DA da;
  Vec a,b,c;
  FFT_DATA *fft_a, *fft_b, *fft_c;
  real x_mid[3], r[3], r_s, h[3];
  int x[3], n[3], i[3], max_i, N[3];
  PetscScalar ***a_vec, ***b_vec, L;
  int index=0;

  PD = params->PD;
  da = params->da;
  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    x_mid[dim] = PD->interval[0]+0.5*(PD->interval[1]-PD->interval[0]);
  FOR_DIM
    N[dim] = PD->N[dim];
  L= PD->interval[1]-PD->interval[0];


  /* get size of local grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  max_i = n[0]*n[1]*n[2];

  VecDuplicate(params->x, &a);
  VecDuplicate(params->x, &b);
  VecDuplicate(params->x, &c);

  DAVecGetArray(da, a, &a_vec);
  DAVecGetArray(da, b, &b_vec);

   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+PD->interval[0]-x_mid[dim];

	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  a_vec[i[2]][i[1]][i[0]] = exp(-r_s*r_s/2.)/pow(2.*M_PI,1.5);

	  FOR_DIM
	    {
	      r[dim] = i[dim]*h[dim];
	      if( i[dim]>=N[dim]/2)
		r[dim] -= L;
	    }
	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  b_vec[i[2]][i[1]][i[0]] = exp(-r_s*r_s/2.)/pow(2.*M_PI,1.5);
	}

  DAVecRestoreArray(da, a, &a_vec);
  DAVecRestoreArray(da, b, &b_vec);

  fft_a=NULL;
  fft_b=NULL;
  fft_a =  ComputeFFTfromVec(da, params->fft_plan, a, fft_a);
  fft_b =  ComputeFFTfromVec(da, params->fft_plan, b, fft_b);
  fft_c = (FFT_DATA*) calloc(max_i, sizeof(*fft_c));

  for(index=0; index<max_i; index++)
    {
      fft_c[index].re = fft_a[index].re*fft_b[index].re
	- fft_a[index].im*fft_b[index].im;
      fft_c[index].im = fft_a[index].re*fft_b[index].im
	+ fft_a[index].im*fft_b[index].re;

    }

  ComputeVecfromFFT(da, params->fft_plan, c, fft_c);
  VecScale(c,  PD->h[0]*PD->h[1]*PD->h[2]/PD->N[0]/PD->N[1]/PD->N[2]);

  VecView(a,PETSC_VIEWER_STDERR_WORLD);
  VecView(b,PETSC_VIEWER_STDERR_WORLD);
  VecView(c,PETSC_VIEWER_STDERR_WORLD);



  VecDestroy(a);
  VecDestroy(b);
  VecDestroy(c);

  free(fft_a);
  free(fft_b);
  free(fft_c);
}


#ifdef MATPRECOND
void TestPreconditioner(MatPrecond MP, Vec x, Vec y)
{
  Mat P;

  P = MP->P;

  KSPSolve(MP->ksp, x, y);

  VecView(x,PETSC_VIEWER_STDERR_WORLD);
  VecView(y,PETSC_VIEWER_STDERR_WORLD);
  MatView(P,PETSC_VIEWER_STDERR_WORLD);
  exit(1);

}
#endif




