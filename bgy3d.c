/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3d.c,v 1.48 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-vec.h"          /* bgy3d_vec_destroy() */
#include "bgy3d-dirichlet.h"    /* bgy3d_laplace_create() */

/* Set  on startup  in bgy3d-main.c.   Used read-only  in a  few other
   files.  Moved here  because  bgy3d-main.c is  not  linked when  the
   main() is taken from elsewhere: */
int verbosity = 0;

/* Get  problem data  (e.g.  from  command line)  using bgy3d_getopt_*
   interface: */
ProblemData bgy3d_problem_data (void)
{
    ProblemData PD;
    int N = 0;

    /* defaults: for older solvers: */
    /* real maxL = 5.0; */
    /* real beta = 0.6061; */

    /* default for recent solvers: */
    real maxL = 12.0;
    real beta = 1.6889;

    real rho = 0.3, h = 0.5;

    /* Grid points in 1 dimension */
    bgy3d_getopt_int ("--N", &N);

    /* inverse temperature */
    bgy3d_getopt_real ("--beta", &beta);

    /* Density */
    bgy3d_getopt_real ("--rho", &rho);

    /* mesh width */
    bgy3d_getopt_real ("--mesh", &h);

    /* (half of the) box size: */
    bgy3d_getopt_real ("--L", &maxL);

    /*=================================*/
    /* set Problem Data */
    /*=================================*/

    /* It appears the intervals for x, y, and z are the same: */
    PD.interval[0] = -maxL;
    PD.interval[1] = +maxL;

    if (N == 0) {
        assert (h > 0.0);
        N = (int) ceil((PD.interval[1] - PD.interval[0]) / h);
    }
    else {
        assert (N > 0);
        h = (PD.interval[1] - PD.interval[0]) / N;
    }
    /* At  this  point  N  and   h  have  consistent  values  for  the
       interval. */

    /* FIXME:  N^2 +  N^2 +  N^2 should  not overflow,  this condition
       ensures that 3N^2 < 2^31: */
    assert (N < 26755);

    FOR_DIM
        PD.N[dim] = N;

    /* FIXME: N^3 should not overflow, this condition ensures that N^3
       < 2^31: */
    assert (N < 1291);

    PD.N3 = PD.N[0] * PD.N[1] * PD.N[2];

    FOR_DIM
        PD.h[dim] = h;

    PD.beta = beta;
    PD.rho  = rho;

    /*
     * Other  parameters  were traditionally  handled  by the  solvers
     * themself each  having potentially its own set  of the defaults.
     * Read commonly used flags from command line:
     */

    /* Mixing parameter: */
    PD.lambda = 0.1;
    bgy3d_getopt_real ("--lambda", &PD.lambda);

    assert (PD.lambda >= 0.0);
    assert (PD.lambda <= 1.0);

    /* Get damp_start from command line*/
    PD.damp = 0.0;
    bgy3d_getopt_real ("--damp-start", &PD.damp);

    /* Number of total iterations */
    PD.max_iter = 100;
    bgy3d_getopt_int ("--max-iter", &PD.max_iter);

    /* Norm_tol for convergence test */
    PD.norm_tol = 1.0e-2;
    bgy3d_getopt_real ("--norm-tol", &PD.norm_tol);

    /*
      Zeropad.  FIXME:  The  code   appears  to  break  when  zpad  >
      L. Regression tests have them  equal, so make it the default to
      have one command line flag fewer:
    */
    PD.zpad = maxL;
    bgy3d_getopt_real ("--zpad", &PD.zpad);

    return PD;
}


State* bgy3d_state_make (const ProblemData *PD)
{
  /* Also initialize all pointers stored in State to NULL: */
  State *BHD = calloc (1, sizeof *BHD);

  BHD->PD = PD;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &BHD->fft_mat, &BHD->da, &BHD->dc);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  BHD->dirichlet_mat = bgy3d_dirichlet_create (BHD->da, BHD->PD);
#endif

#ifdef L_BOUNDARY_MG
  /* multigrid, apparently needs two descriptors: */
#error "Need BHD->da_dmmg"
#endif

  /* Create global scratch vectors: */
  bgy3d_vec_create1 (BHD->da, 3, BHD->v); /* real */

  /* Complex  vectors for  k-space representations.   These  three are
     used by ComputeFFTfromCoulomb(): */
  bgy3d_vec_create1 (BHD->dc, 3, BHD->fg2_fft); /* complex */

  BHD->fft_scratch = bgy3d_vec_create (BHD->dc); /* complex */

  return BHD;
}

void bgy3d_state_destroy (State *BHD)
{
  MPI_Barrier (PETSC_COMM_WORLD);

  bgy3d_vec_destroy1 (3, BHD->v);
  bgy3d_vec_destroy1 (3, BHD->fg2_fft);

  bgy3d_vec_destroy (&BHD->fft_scratch);

#ifdef L_BOUNDARY
  assert (BHD->dirichlet_mat != NULL);
  MatDestroy (BHD->dirichlet_mat);
#endif

#ifdef L_BOUNDARY_MG
  DMMGDestroy (BHD->dmmg);
#endif

  DADestroy (BHD->da);
  DADestroy (BHD->dc);
  MatDestroy (BHD->fft_mat);

  free (BHD);
}

/* Code used to be verbose: */
void bgy3d_state_print (const State *BHD)
{
  const ProblemData *PD = BHD->PD;
  const real L = PD->interval[1] - PD->interval[0];

  PetscPrintf (PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  PetscPrintf (PETSC_COMM_WORLD, "N = %d x %d x %d\n", PD->N[0], PD->N[1], PD->N[2]);
  PetscPrintf (PETSC_COMM_WORLD, "h = %g x %g x %g\n", PD->h[0], PD->h[1], PD->h[2]);
  PetscPrintf (PETSC_COMM_WORLD, "β = %g (%5.1f K)\n", PD->beta, 1.0 / PD->beta / KBOLTZMANN);
  PetscPrintf (PETSC_COMM_WORLD, "ρ = %g (%g per cell)\n", PD->rho, PD->rho * L * L * L);
}

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

  /* First entry in the file appears to be the number of sites: */
  if (1 != fscanf (fp, "%d", N))
    {
      PetscPrintf (PETSC_COMM_WORLD, "Error reading file (molecule).\n");
      exit (1);
    }

  if(verbosity>2)
    PetscPrintf(PETSC_COMM_WORLD,"Reading molecule data:\n%d particles\n", *N);

  x_M= (real**) malloc(sizeof(*x_M)*(*N));


  for(i=0; i<*N; i++)
    {
      x_M[i] = (real*) malloc(sizeof(real)*3);
      if (3 != fscanf (fp, "%lf %lf %lf\n", &x_M[i][0], &x_M[i][1], &x_M[i][2]))
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

/* Reduce buffer by summing respective entries on all workeres: */
void bgy3d_comm_allreduce (int n, real x[n])
{
  int err = MPI_Allreduce (MPI_IN_PLACE, x, n, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  assert (!err);
}

