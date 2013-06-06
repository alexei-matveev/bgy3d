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


/*
  Get  problem data  (e.g.   from command  line) using  bgy3d_getopt_*
  interface.

  FIXME:  some parameters are  still handled  by the  solvers themself
  each having potentially its own set of the defaults.

  Do not forget to set the  defaults, the command line may not contain
  a setting for every ProblemData field!
*/
ProblemData bgy3d_problem_data (void)
{
    ProblemData PD;

    /* Grid points in 1 dimension */
    int N = 32;
    bgy3d_getopt_int ("--N", &N);
    assert (N > 0);

    /* (half of the) box size: */
    real maxL = 12.0;
    bgy3d_getopt_real ("--L", &maxL);

    /* It appears the intervals for x, y, and z are the same: */
    FOR_DIM
      PD.L[dim] = 2 * maxL;

    FOR_DIM
      PD.N[dim] = N;

    FOR_DIM
      PD.h[dim] = PD.L[dim] / PD.N[dim];

    /*
      At this  point N, h and  L have consistent values.  FIXME: N^2 +
      N^2 + N^2 should not  overflow, this condition ensures that 3N^2
      < 2^31:
    */
    assert (N < 26755);

    /* FIXME: N^3 should not overflow, this condition ensures that N^3
       < 2^31: */
    if (!bgy3d_getopt_test ("--rism"))
        assert (N < 1291);

    /* Inverse temperature: */
    PD.beta = 1.6889;
    bgy3d_getopt_real ("--beta", &PD.beta);

    /* Density: */
    PD.rho = 0.3;
    bgy3d_getopt_real ("--rho", &PD.rho);

    /* Mixing parameter: */
    PD.lambda = 0.1;
    bgy3d_getopt_real ("--lambda", &PD.lambda);

    assert (PD.lambda >= 0.0);
    assert (PD.lambda <= 1.0);

    /*
      Get initial  interaction scaling  factor from command  line. The
      code  used   to  automate  "blending-in"   some  interaction  by
      iterating  scaling factor  from  some low  initial  value to  1.
      Default was  0.0 originally, but that makes  necessary to always
      specify --damp-start when 1.0 you do not want blending.  I think
      1.0 is a more reasonable default:
    */
    PD.damp = 1.0;
    bgy3d_getopt_real ("--damp-start", &PD.damp);

    /* Number of total iterations: */
    PD.max_iter = 100;
    bgy3d_getopt_int ("--max-iter", &PD.max_iter);

    /* Norm tolerance for convergence test: */
    PD.norm_tol = 1.0e-2;
    bgy3d_getopt_real ("--norm-tol", &PD.norm_tol);

    /*
      Zeropad.  FIXME:  The  code   appears  to  break  when  zpad  >
      L. Regression tests have them  equal, so make it the default to
      have one command line flag fewer:
    */
    PD.zpad = maxL;
    bgy3d_getopt_real ("--zpad", &PD.zpad);


    /* Closure is  only used in  HNC-like methods. Supply  default for
       the case the command line does not provide it: */
    char closure[] = "HNC";
    bgy3d_getopt_string ("--closure", closure, sizeof closure);

    if (strcmp (closure, "HNC") == 0)
      PD.closure = CLOSURE_HNC;
    else if (strcmp (closure, "KH") == 0)
      PD.closure = CLOSURE_KH;
    else if (strcmp (closure, "PY") == 0)
      PD.closure = CLOSURE_PY;
    else
      {
        PetscPrintf (PETSC_COMM_WORLD, "No such OZ closure: %s\n", closure);
        exit (1);
      }

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

  return BHD;
}

void bgy3d_state_destroy (State *BHD)
{
  MPI_Barrier (PETSC_COMM_WORLD);

#ifdef L_BOUNDARY
  assert (BHD->dirichlet_mat != NULL);
  MatDestroy (&BHD->dirichlet_mat);
#endif

#ifdef L_BOUNDARY_MG
  DMMGDestroy (BHD->dmmg);
#endif

  DMDestroy (&BHD->da);
  DMDestroy (&BHD->dc);
  MatDestroy (&BHD->fft_mat);

  free (BHD);
}

/* Code used to be verbose: */
void bgy3d_problem_data_print (const ProblemData *PD)
{
  const real L3 = PD->L[0] * PD->L[1] * PD->L[2];

  PetscPrintf (PETSC_COMM_WORLD, "Ω = %g x %g x %g A³\n", PD->L[0], PD->L[1], PD->L[2]);
  PetscPrintf (PETSC_COMM_WORLD, "N = %d x %d x %d\n", PD->N[0], PD->N[1], PD->N[2]);
  PetscPrintf (PETSC_COMM_WORLD, "h = %g x %g x %g A³\n", PD->h[0], PD->h[1], PD->h[2]);
  PetscPrintf (PETSC_COMM_WORLD, "β = %g (%5.1f K)\n", PD->beta, 1.0 / PD->beta / KBOLTZMANN);
  PetscPrintf (PETSC_COMM_WORLD, "ρ = %g (%g per cell)\n", PD->rho, PD->rho * L3);
  PetscPrintf (PETSC_COMM_WORLD, "a = %g (Seitz radius)\n", pow ((4 * M_PI / 3) * PD->rho, -1.0 / 3.0));
  PetscPrintf (PETSC_COMM_WORLD, "λ = %g (mixing ratio)\n", PD->lambda);
  PetscPrintf (PETSC_COMM_WORLD, "norm-tol = %e\n", PD->norm_tol);
  PetscPrintf (PETSC_COMM_WORLD, "max-iter = %d\n", PD->max_iter);
  PetscPrintf (PETSC_COMM_WORLD, "zpad = %g\n", PD->zpad);
  switch (PD->closure)
    {
    case CLOSURE_HNC:
      PetscPrintf (PETSC_COMM_WORLD, "closure = HNC\n");
      break;
    case CLOSURE_KH:
      PetscPrintf (PETSC_COMM_WORLD, "closure = KH\n");
      break;
    case CLOSURE_PY:
      PetscPrintf (PETSC_COMM_WORLD, "closure = PY\n");
      break;
    }
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

