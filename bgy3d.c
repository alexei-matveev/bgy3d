/*==========================================================*/
/*  $Id: bgy3d.c,v 1.48 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"

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
