/*==========================================================*/
/*  $Id: bgy3d.c,v 1.48 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"

/* Get  problem data  (e.g.  from  command line)  using bgy3d_getopt_*
   interface: */
ProblemData bgy3d_problem_data (void)
{
    ProblemData PD;
    int N = 0;
    real beta = 0.6061, rho = 0.3, h = 0.5, interval[2] = {-5.0, 5.0};

    /* Grid points in 1 dimension */
    bgy3d_getopt_int ("-N", &N);

    /* inverse temperature */
    bgy3d_getopt_real ("-beta", &beta);

    /* Density */
    bgy3d_getopt_real ("-rho", &rho);

    /* mesh width */
    bgy3d_getopt_real ("-mesh", &h);

    /*=================================*/
    /* set Problem Data */
    /*=================================*/

    if (N == 0)
        N = (int) ceil((interval[1] - interval[0]) / h);
    else
        h = (interval[1] - interval[0]) / N;

    FOR_DIM
        PD.N[dim] = N;

    PD.N3 = PD.N[0] * PD.N[1] * PD.N[2];

    FOR_DIM
        PD.h[dim] = h;

    PD.interval[0] = interval[0];
    PD.interval[1] = interval[1];

    PD.beta = beta;
    PD.rho  = rho;

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
