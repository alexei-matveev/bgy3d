/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3d.c,v 1.48 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"

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

    FOR_DIM
        PD.N[dim] = N;

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

    assert (PD.lambda > 0.0);
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

    /* Parallel staff: */
    MPI_Comm_size(PETSC_COMM_WORLD, &PD.np);
    MPI_Comm_rank(PETSC_COMM_WORLD, &PD.id);

    return PD;
}

real Lennard_Jones(real r, real epsilon, real sigma)
{
  real sr6, sr, re;

  r=r+SHIFT;

  sr = sigma/r;
  sr6 = SQR(sr)*SQR(sr)*SQR(sr);

  re= 4.*epsilon*sr6*(sr6-1.);

  if(fabs(re)>epsilon*CUTOFF)
    return epsilon*CUTOFF;
  else
    return re;
}

real Lennard_Jones_grad(real r, real xr, real epsilon, real sigma)
{
  real sr6, sr, re;

  r=r+SHIFT;

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
