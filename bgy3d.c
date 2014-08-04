/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-vec.h"          /* bgy3d_vec_destroy() */
#include "bgy3d-dirichlet.h"    /* bgy3d_laplace_create() */

/*
  Set  on  startup  by  bgy3d_guile_init()  in  bgy3d-guile.c  and  in
  bgy3d-main.c.   Used read-only  in a  few other  files.   Moved here
  because bgy3d-main.c  is not  linked when the  main() is  taken from
  elsewhere:
*/
int verbosity = 0;

/* This  communicator  will be  flipped  between PETSC_COMM_WORLD  and
   PETSC_COMM_SELF occasionally: */
MPI_Comm comm_world_petsc = MPI_COMM_NULL;


State* bgy3d_state_make (const ProblemData *PD)
{
  /*
    FIXME: Memory limits? Accidentally  calling 3D code with a typical
    dimension  of   1D?   Anyway,  N^3  should  not   overflow  in  3D
    runs. Print a warning if N^3 is definitely over 2^31:
  */
  const int nmax = 1291;
  if (PD->N[0] >= nmax && PD->N[1] >= nmax && PD->N[2] >= nmax)
    fprintf (stderr, "Warning: grid %d x %d x %d too large for 3D!\n",
             PD->N[0], PD->N[1], PD->N[2]);

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
  MPI_Barrier (comm_world_petsc);

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
  const real L3 = volume (PD);

  PRINTF ("Ω = %g x %g x %g A³\n", PD->L[0], PD->L[1], PD->L[2]);
  PRINTF ("N = %d x %d x %d\n", PD->N[0], PD->N[1], PD->N[2]);
  PRINTF ("h = %g x %g x %g A³\n", PD->h[0], PD->h[1], PD->h[2]);
  PRINTF ("β = %g kcal⁻¹ (%5.1f K)\n", PD->beta, 1.0 / PD->beta / KBOLTZMANN);
  PRINTF ("ρ = %g A⁻³ (%g per cell)\n", PD->rho, PD->rho * L3);
  PRINTF ("a = %g A (Seitz radius)\n", pow ((4 * M_PI / 3) * PD->rho, -1.0 / 3.0));
  PRINTF ("λ = %g (mixing ratio)\n", PD->lambda);
  PRINTF ("norm-tol = %e\n", PD->norm_tol);
  PRINTF ("max-iter = %d\n", PD->max_iter);

  switch (PD->closure)
    {
    case CLOSURE_HNC:
      PRINTF ("closure = HNC\n");
      break;
    case CLOSURE_KH:
      PRINTF ("closure = KH\n");
      break;
    case CLOSURE_PY:
      PRINTF ("closure = PY\n");
      break;
    case CLOSURE_PSE0:
    case CLOSURE_PSE1:
    case CLOSURE_PSE2:
    case CLOSURE_PSE3:
    case CLOSURE_PSE4:
    case CLOSURE_PSE5:
    case CLOSURE_PSE6:
    case CLOSURE_PSE7:
      PRINTF ("closure = PSE%d\n", PD->closure - CLOSURE_PSE0);
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
      PRINTF ("File not found.\n");
      exit(1);
    }

  /* First entry in the file appears to be the number of sites: */
  if (1 != fscanf (fp, "%d", N))
    {
      PRINTF ("Error reading file (molecule).\n");
      exit (1);
    }

  if(verbosity>2)
    PRINTF ("Reading molecule data:\n%d particles\n", *N);

  x_M= (real**) malloc(sizeof(*x_M)*(*N));


  for(i=0; i<*N; i++)
    {
      x_M[i] = (real*) malloc(sizeof(real)*3);
      if (3 != fscanf (fp, "%lf %lf %lf\n", &x_M[i][0], &x_M[i][1], &x_M[i][2]))
	{
	  PRINTF ("Error reading file (molecule).\n");
	  exit(1);
	}
      if(verbosity>2)
	PRINTF ("%f %f %f\n", x_M[i][0],
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


/*******************/
/* COMM PRIMITIVES */
/*******************/

/* Reduce buffer by summing respective entries on all workeres: */
void comm_allreduce (int n, real x[n])
{
  int err =
    MPI_Allreduce (MPI_IN_PLACE, x, n, MPI_DOUBLE, MPI_SUM, comm_world_petsc);
  assert (!err);
}


/* Return MPI runk in comm_world_petsc: */
int comm_rank (void)
{
  int rank;
  MPI_Comm_rank (comm_world_petsc, &rank);
  return rank;
}


/* Return MPI size of comm_world_petsc: */
int comm_size (void)
{
  int size;
  MPI_Comm_size (comm_world_petsc, &size);
  return size;
}


/*
  If the  flag is #f set  the parallel mode, otherwise  set the serial
  mode of  operation for PETSC. The  current mode is  returned and may
  have to be restored later.
*/
bool
comm_set_parallel_x (bool flag)
{
  const MPI_Comm comm = comm_world_petsc;
  assert (comm == PETSC_COMM_WORLD || comm == PETSC_COMM_SELF);

  MPI_Barrier (PETSC_COMM_WORLD);
  if (flag)
    comm_world_petsc = PETSC_COMM_WORLD;
  else
    comm_world_petsc = PETSC_COMM_SELF;

  return (comm == PETSC_COMM_WORLD);
}
