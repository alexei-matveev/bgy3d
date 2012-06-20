/*==========================================================*/
/*  $Id: bgy3dH2O.c,v 1.42 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"


#include "bgy3d_SolventParameters.h"




#define R_r  9
#define R_l  7
#define ITERTOL_DOWN 0.9
#define ITERTOL_UP   5.0e-1


#define ZEROPAD 5

#define SL 16.0
#define SR 16.0
//#define NORM_REG 1.0e-1
//#define NORM_REG 8.0e-1

real  NORM_REG=1.0e-1;
real  NORM_REG2=1.0e-2;

BGY3dH2OData BGY3dH2OData_Pair_malloc(PData PD)
{
  BGY3dH2OData BHD;
  DA da;
  real interval[2], h[3], N[3], L, beta, maxL;
  int x[3], n[3], dim;

  int np;
  int local_nx, local_x_start, local_ny, local_y_start, total_local_size;
  PetscInt lx[1], ly[1], *lz;


  BHD = (BGY3dH2OData) malloc(sizeof(*BHD));

  /****************************************************/
  /* set Lennard-Jones and Coulomb parameters */
  /****************************************************/

  /* water hydrogen */
  // BHD->LJ_paramsH = (void* ) malloc(sizeof(real)*3);
  // ((real*)(BHD->LJ_paramsH))[0] = eH;   /* espilon */
  // ((real*)(BHD->LJ_paramsH))[1] = sH;   /* sigma   */
  // ((real*)(BHD->LJ_paramsH))[2] = SQR(qH);   /* q   */
  BHD->LJ_paramsH[0] = eH;  /* epsilon  */
  BHD->LJ_paramsH[1] = sH;  /* sigma    */
  BHD->LJ_paramsH[2] = SQR(qH); /* charge product */

  /* water oxygen */
  // BHD->LJ_paramsO = (void* ) malloc(sizeof(real)*3);
  // ((real*)(BHD->LJ_paramsO))[0] = eO;   /* espilon */
  // ((real*)(BHD->LJ_paramsO))[1] = sO;   /* sigma   */
  // ((real*)(BHD->LJ_paramsO))[2] = SQR(qO);   /* q   */
  BHD->LJ_paramsO[0] = eO;  /* epsilon  */
  BHD->LJ_paramsO[1] = sO;  /* sigma    */
  BHD->LJ_paramsO[2] = SQR(qO); /* charge product */

  /* water O-H mixed parameters */
  // BHD->LJ_paramsHO = (void* ) malloc(sizeof(real)*3);
  // ((real*)(BHD->LJ_paramsHO))[0] = sqrt(eH*eO);   /* espilon */
  // ((real*)(BHD->LJ_paramsHO))[1] = 0.5*(sH+sO);  /* sigma   */
  // ((real*)(BHD->LJ_paramsHO))[2] = qH*qO;         /* q   */
  BHD->LJ_paramsHO[0] = sqrt(eH*eO);  /* epsilon  */
  BHD->LJ_paramsHO[1] = 0.5*(sH+sO);  /* sigma    */
  BHD->LJ_paramsHO[2] = qH*qO; /* charge product */

  /****************************************************/


  BHD->PD = PD;
  /*****************************/
  /* reset standard parameters */
  /*****************************/
  maxL=12.0;
  PetscOptionsGetReal(PETSC_NULL,"-L",&maxL, PETSC_NULL);
  PD->interval[0] = -maxL;//-25.0;
  PD->interval[1] = maxL;//25.0;
  FOR_DIM
    PD->h[dim] = (PD->interval[1]-PD->interval[0])/PD->N[dim];
  PD->N3 = PD->N[0]*PD->N[1]*PD->N[2];
  beta= 1.6889;
  PetscOptionsGetReal(PETSC_NULL,"-beta",&beta, PETSC_NULL);
  PD->beta = beta;//1.6889;
  PetscPrintf(PETSC_COMM_WORLD, "Corrected domain size:\n");
  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  //PetscPrintf(PETSC_COMM_WORLD, "Boundary smoothing parameters : SL= %f  SR= %f\n", SL, SR);
  //PetscPrintf(PETSC_COMM_WORLD, "ZEROPAD= %f\n", ZEROPAD);
  PetscPrintf(PETSC_COMM_WORLD, "Regularization of normalization: NORM_REG= %e\n", NORM_REG);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", PD->beta);
  /******************************/
  PetscPrintf(PETSC_COMM_WORLD,"--------------------------------\n");
  PetscPrintf(PETSC_COMM_WORLD,"eps1= %f \t eps2= %f\n", eH, eO);
  PetscPrintf(PETSC_COMM_WORLD,"sig1= %f \t sig2= %f\n", sH, sO);
  PetscPrintf(PETSC_COMM_WORLD,"q1  = %f \t q2  = %f\n", qH, qO);
  PetscPrintf(PETSC_COMM_WORLD,"--------------------------------\n");

  BHD->beta = PD->beta;
  BHD->rho  = PD->rho;
  BHD->rho_H = PD->rho; //2.*PD->rho;
  BHD->rho_O = PD->rho;
  beta = PD->beta;

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];
  L=interval[1]-interval[0];
  FOR_DIM
    h[dim]=PD->h[dim];
  FOR_DIM
    N[dim]=PD->N[dim];

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
  DACreateGlobalVector(da, &(BHD->xHO));
  VecSet(BHD->xH, 0.0);
  VecSet(BHD->xO, 0.0);
  VecSet(BHD->xHO, 0.0);
#else
  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR ,
	     PD->N[0], PD->N[1], PD->N[2],
	     1, 1, np,
	     1,0,
	     lx, ly, lz,
	     &(BHD->da));


  da = BHD->da;
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
  /* Create global vectors */
  DACreateGlobalVector(da, &(BHD->gH_ini));
  DACreateGlobalVector(da, &(BHD->gO_ini));
  DACreateGlobalVector(da, &(BHD->gHO_ini));
  DACreateGlobalVector(da, &(BHD->ucH));
  DACreateGlobalVector(da, &(BHD->ucO));
  DACreateGlobalVector(da, &(BHD->ucHO));
  DACreateGlobalVector(da, &(BHD->cH));
  DACreateGlobalVector(da, &(BHD->cO));
  DACreateGlobalVector(da, &(BHD->cHO));
  FOR_DIM
    {
/*       VecDuplicate(BHD->gH_ini, &(BHD->fH[dim])); */
/*       VecDuplicate(BHD->gH_ini, &(BHD->fO[dim])); */
/*       VecDuplicate(BHD->gH_ini, &(BHD->fHO[dim])); */
/*       VecDuplicate(BHD->gH_ini, &(BHD->fH_l[dim])); */
/*       VecDuplicate(BHD->gH_ini, &(BHD->fO_l[dim])); */
/*       VecDuplicate(BHD->gH_ini, &(BHD->fHO_l[dim])); */
/*       VecDuplicate(BHD->gH_ini, &(BHD->v[dim])); */

      DACreateGlobalVector(da, &(BHD->fH[dim]));
      DACreateGlobalVector(da, &(BHD->fO[dim]));
      DACreateGlobalVector(da, &(BHD->fHO[dim]));
      DACreateGlobalVector(da, &(BHD->fH_l[dim]));
      DACreateGlobalVector(da, &(BHD->fO_l[dim]));
      DACreateGlobalVector(da, &(BHD->fHO_l[dim]));
      DACreateGlobalVector(da, &(BHD->v[dim]));
    }




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
      BHD->fg2_fft[dim] = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
    }

  BHD->g_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->gfg2_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->fft_scratch = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->ucH_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->ucO_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->ucHO_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
  BHD->wHO_fft = (fftw_complex*) malloc(n[0]*n[1]*n[2]*sizeof(fftw_complex));
   /* Compute initial data */
  RecomputeInitialData(BHD, 1.0, 1.0);

  free(lz);

  return BHD;
}



void BGY3dH2OData_free(BGY3dH2OData BHD)
{
  int dim;

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
      free(BHD->fg2_fft[dim]);
    }
  free(BHD->g_fft);
  free(BHD->gfg2_fft);
  free(BHD->fft_scratch);
  free(BHD->ucH_fft);
  free(BHD->ucO_fft);
  free(BHD->ucHO_fft);
  free(BHD->wHO_fft);

  VecDestroy(BHD->gH_ini);
  VecDestroy(BHD->gO_ini);
  VecDestroy(BHD->gHO_ini);
  VecDestroy(BHD->ucH);
  VecDestroy(BHD->ucO);
  VecDestroy(BHD->ucHO);
  VecDestroy(BHD->cH);
  VecDestroy(BHD->cO);
  VecDestroy(BHD->cHO);
  DADestroy(BHD->da);
#ifdef L_BOUNDARY
  MatDestroy(BHD->M);
  KSPDestroy(BHD->ksp);
  VecDestroy(BHD->xH);
  VecDestroy(BHD->xHO);
  VecDestroy(BHD->xO);
#endif
  // No need to do this anymore
  // free(BHD->LJ_paramsH);
  // free(BHD->LJ_paramsO);
  // free(BHD->LJ_paramsHO);

  fftwnd_mpi_destroy_plan(BHD->fft_plan_fw);
  fftwnd_mpi_destroy_plan(BHD->fft_plan_bw);

  free(BHD);
}

// real LJ_repulsive(real r, void *LJ_params)
real LJ_repulsive(real r, real epsilon, real sigma)
{
  // real epsilon, sigma, sr6, sr, re;
  real sr6, sr, re;

  r=r+SHIFT;

  // epsilon = ((double*)LJ_params)[0];
  // sigma   = ((double*)LJ_params)[1];

  sr = sigma/r;
  sr6 = SQR(sr)*SQR(sr)*SQR(sr);

  re= 4.*epsilon*sr6*sr6;

  if(fabs(re)>epsilon*CUTOFF)
    return epsilon*CUTOFF;
  else
    return re;
}


// Alternate void pointer as real number
// real Coulomb_short( real r, void *params)
real Coulomb_short( real r, real q2)
{
  // real q2, re;
  real re;

  // q2 = ((double*)params)[2];

   if(r==0)
     {
       if( q2>0)
	 return EPSILON0INV * q2 * (CUTOFF*1.0e-5);
       else
	 return EPSILON0INV * q2 * (CUTOFF*1.0e-5); //1.0e+4;
     }



  re = EPSILON0INV * q2 * erfc(G*r)/r;



  /* check for large values */
  /* remember: exp(-re) will be computed */
  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return EPSILON0INV * q2 * (CUTOFF*1.0e-5);
/*   else if( -re>1e+1) */
/*     return -1e+1; */
  else
    return re;
}

// Alternate void pointer as real number
// real Coulomb_short_grad( real r, real rx, void *params)
real Coulomb_short_grad( real r, real rx, real q2 )
{
  // real q2, re;
  real re;

  // q2 = ((double*)params)[2];

  if(rx==0)
    return 0;
  if(r==0)
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);


  re = - EPSILON0INV * q2 * (erfc(G*r) + 2.*G/sqrt(M_PI)*r*exp(-G*G*r*r))*rx/pow(r,3.0);

  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);
  else
    return re;
}

// real Coulomb_long( real r, void *params)
real Coulomb_long( real r, real q2)
{
  // real q2, re;
  real re;

  // q2 = ((double*)params)[2];

   if(r==0)
     {
       return EPSILON0INV * q2 * G * 2./sqrt(M_PI);
     }

  re = EPSILON0INV * q2 * erf(G*r)/r;



  /* check for large values */
  /* remember: exp(-re) will be computed */
  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return EPSILON0INV * q2 * (CUTOFF*1.0e-5);
/*   else if( -re>1e+1) */
/*     return -1e+1; */
  else
    return re;
}

// real Coulomb_long_grad( real r, real rx, void *params)
real Coulomb_long_grad( real r, real rx, real q2)
{
  // real q2, re;
  real re;

  // q2 = ((double*)params)[2];

  if(r==0)
    return 0;


  re = - EPSILON0INV * q2 * (erf(G*r)
			     - 2.*G/sqrt(M_PI)*r*exp(-G*G*r*r))*rx/pow(r,3.0);

  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);
  else
    return re;
}

// real Coulomb_long_spline( real r, void *params)
real Coulomb_long_spline( real r, real q2)
{
  // real q2, re;
  real re;
  real r_rl_2, rr_rl_2r, rr_rl_3, s;

  // q2 = ((double*)params)[2];

  if(r==0)
    {
      return EPSILON0INV * q2 * G * 2./sqrt(M_PI);
    }

  if( r > R_r)
    re = 0;
  else if ( r > R_l)
    {
      r_rl_2 = SQR(r-R_l);
      rr_rl_2r = (3.*R_r-R_l-2.*r);
      rr_rl_3 = pow(R_r-R_l,3);

      s= 1. - r_rl_2*rr_rl_2r/rr_rl_3;

      re = s * EPSILON0INV * q2 * erf(G*r)/r;

    }
  else
    re = EPSILON0INV * q2 * erf(G*r)/r;



  /* check for large values */
  /* remember: exp(-re) will be computed */
  if(fabs(re) > fabs(EPSILON0INV * q2 * CUTOFF))
    return EPSILON0INV * q2 * CUTOFF;
/*   else if( -re>1e+1) */
/*     return -1e+1; */
  else
    return re;
}

// real Coulomb_long_spline_grad( real r, real rx, void *params)
real Coulomb_long_spline_grad( real r, real rx, real q2)
{
  // real q2, re;
  real re;
  real r_rl, rr_rl_2r, rr_rl_3, s, ss;


  // q2 = ((double*)params)[2];

  if(r==0)
    return 0;

  if( r > R_r)
    re = 0;
  else if ( r > R_l)
    {
      r_rl = (r-R_l);
      rr_rl_2r = (3.*R_r-R_l-2.*r);
      rr_rl_3 = pow(R_r-R_l,3);

      s= 1. - SQR(r_rl)*rr_rl_2r/rr_rl_3;
      ss = 2.*r_rl*(r_rl - rr_rl_2r)/rr_rl_3;

      re = - EPSILON0INV * q2 *
	(s * ((erf(G*r) - 2.*G/sqrt(M_PI)*r*exp(-G*G*r*r))*rx/pow(r,3.0))
	 - ss * erf(G*r) * rx/SQR(r)) ;

    }
  else
    re = -  EPSILON0INV * q2 * (erf(G*r)
			       - 2.*G/sqrt(M_PI)*r*exp(-G*G*r*r))*rx/pow(r,3.0);


  if(fabs(re) > fabs(EPSILON0INV * q2 * CUTOFF))
    return -EPSILON0INV * q2 * CUTOFF;
  else
    return re;
}


// real Coulomb( real r, void *params)
real Coulomb( real r, real q2)
{
  // real q2, re;
  real  re;

  // q2 = ((double*)params)[2];

   if(r==0)
     {

       return EPSILON0INV * q2 * (CUTOFF*1.0e-5);


     }


  re = EPSILON0INV * q2 /r;

  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)))
    return EPSILON0INV * q2 * (CUTOFF*1.0e-5);
/*   else if( -re>1e+1) */
/*     return -1e+1; */
  else
    return re;
}

// real Coulomb_grad( real r, real rx, void *params)
real Coulomb_grad( real r, real rx, real q2)
{
  // real q2, re;
  real re;

  // q2 = ((double*)params)[2];

  if( rx==0 )
    return 0;
  if(r==0)
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);

  re = - EPSILON0INV * q2 * rx/pow(r,3.0);



  if(fabs(re) > fabs(EPSILON0INV * q2 * (CUTOFF*1.0e-5)) )
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);
  else
    return re;
}


void ComputeH2O_g(Vec g, Vec g0, Vec dg)
{
  int local_size, i;
  PetscScalar *g_vec, *dg_vec, *g0_vec;
  real  k; // g_norm



  VecGetArray( g, &g_vec);
  VecGetArray( g0, &g0_vec);
  VecGetArray( dg, &dg_vec);
  VecGetLocalSize(g, &local_size);

  for(i=0; i<local_size; i++)
    {

      k = -g0_vec[i]-dg_vec[i];
      g_vec[i] = exp(k);

      assert(!isinf(g_vec[i]) && !isnan(g_vec[i]));
    }

  VecRestoreArray(g, &g_vec);
  VecRestoreArray(g, &g0_vec);
  VecRestoreArray(dg, &dg_vec);

/*   VecView(g,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */
}

void CheckMax( Vec g, char name[5], real max)
{
  PetscScalar *g_vec;
  real k;
  int local_size, i, counter=0;


  VecGetArray( g, &g_vec);
  VecGetLocalSize(g, &local_size);

  for(i=0; i<local_size; i++)
    {
      k = g_vec[i];
      if( k>max)
	{
	  counter++;
	  g_vec[i]=max;
	}
    }
  VecRestoreArray(g, &g_vec);
  if( counter>0)
    PetscPrintf(PETSC_COMM_WORLD,"Corrected %s %d times.\n", name, counter);

}


void VecSetRandom_H2O(Vec g, real mag)
{
  int local_size, i;
  PetscScalar *g_vec;

  VecGetArray( g, &g_vec);
  VecGetLocalSize(g, &local_size);

  for(i=0; i<local_size; i++)
    {
      g_vec[i] = mag*2.*(0.5-rand()/(real)RAND_MAX);
    }
  VecRestoreArray(g, &g_vec);
}

void ImposeBoundaryCondition_Initialize( BGY3dH2OData BHD, real zpad)
{
  DA da;
  PData PD;
  int x[3], n[3], dim, index, p_id, p_idr=0, *recv_count;


  PD = BHD->PD;
  da = BHD->da;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  index = (int) ceil((-PD->interval[0]-zpad)/PD->h[0]);
  recv_count = (int*) malloc(PD->np*sizeof(int));
  for(dim=0; dim<PD->np; dim++)
    recv_count[dim] = 1;

  if( index      >= x[0] && index      < x[0]+n[0] &&
      PD->N[1]/2 >= x[1] && PD->N[1]/2 < x[1]+n[1] &&
      PD->N[2]/2 >= x[2] && PD->N[2]/2 < x[2]+n[2]   )
    p_id = PD->id;
  else
    p_id = 0;

  MPI_Reduce( &p_id, &p_idr, 1,MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD );
  MPI_Bcast ( &p_idr, 1, MPI_INT, 0, PETSC_COMM_WORLD );


  BHD->p_id = p_idr;
  BHD->p_index = index;
  PetscPrintf(PETSC_COMM_WORLD,"Root id is %d (%d,%d,%d).\n",
	      p_idr, index, PD->N[1]/2, PD->N[2]/2);
  //intf("%d Root id is %d (%d,%d,%d).\n", PD->id, p_idr, index, PD->N[1]/2, PD->N[2]/2);

  free(recv_count);
}


void ImposeBoundaryCondition( BGY3dH2OData BHD, Vec g)
{
  DA da;
  PData PD;
  int x[3], n[3], index;
  real b_data;
  PetscScalar ***g_vec;

  PD = BHD->PD;
  da = BHD->da;

  index = BHD->p_index;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if( PD->id  == BHD->p_id )
    {
      //printf("%d ich  bins %d", PD->id , BHD->p_id);fflush(stdout);
      DAVecGetArray(da, g, &g_vec);
      b_data = g_vec[PD->N[2]/2][PD->N[1]/2][index];
      DAVecRestoreArray(da, g, &g_vec);

    }

  MPI_Bcast ( &b_data, 1, MPI_INT, BHD->p_id, PETSC_COMM_WORLD );

  VecShift(g, -b_data);
}


real ImposeBoundaryConditionII( BGY3dH2OData BHD, Vec g, real zpad)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], index=0, index_sum, dim;
  real b_data, r[3], r_s, h[3], interval[2], data_sum;
  PetscScalar ***g_vec;

  PD = BHD->PD;
  da = BHD->da;

  interval[0]=PD->interval[0];
  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  DAVecGetArray(da, g, &g_vec);
  b_data=0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0];

	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  if( fabs(r_s-zpad)<=h[0]/2)
	    {
	      b_data += g_vec[i[2]][i[1]][i[0]];
	      index++;
	    }
	}
  DAVecRestoreArray(da, g, &g_vec);

  MPI_Reduce( &b_data, &data_sum, 1,MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD );
  MPI_Reduce( &index, &index_sum, 1,MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD );

  if( PD->id==0)
    data_sum /= (real)index_sum;


  MPI_Bcast ( &data_sum, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD );

  //VecShift(g, -b_data);
  return data_sum;
}


void ComputeH2O_Renormalization( BGY3dH2OData BHD, Vec g)
{
  PData PD;
  real vsum, h3, *h;
  // PetscScalar *g_vec;

  PD = BHD->PD;
  h = PD->h;

  h3 = h[0]*h[1]*h[2]/pow(PD->interval[1]-PD->interval[0],3);

  VecSum(g, &vsum);
  //fprintf(stderr,"%e ", vsum*h3);
  VecScale( g, 1./(h3*vsum));

/*   VecGetArray( g, &g_vec); */
/*   fprintf(stderr,"%e ", g_vec[0]); */
/*   VecShift(g,-g_vec[0]); */
/*   VecRestoreArray(g, &g_vec); */
}

void Smooth_Function(BGY3dH2OData BHD, Vec g, real RL, real RR, real shift)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], dim;
  PetscScalar ***g_vec;
  real r[3], r_s, h[3], interval[2], s, r_rl_2, rr_rl_2r, rr_rl_3;

  PD = BHD->PD;
  da = BHD->da;

  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, g, &g_vec);
   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

/* 	  if(i[2]==0 || i[2]==n[2]-1 || i[1]==0 || i[1]==n[1]-1 ||  */
/* 	     i[0]==0 || i[0]==n[0]-1) */
/* 	    g_vec[i[2]][i[1]][i[0]]=0; */

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0];
	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

	  if( r_s > RR)
	    s=0;
	  else if ( r_s > RL)
	    {
	      r_rl_2 = SQR(r_s-RL);
	      rr_rl_2r = (3.*RR-RL-2.*r_s);
	      rr_rl_3 = pow(RR-RL,3);

	      s= 1. - r_rl_2*rr_rl_2r/rr_rl_3;

	    }
	  else
	    s=1;
	  g_vec[i[2]][i[1]][i[0]] = s* (g_vec[i[2]][i[1]][i[0]]-shift) + shift;
	}
  DAVecRestoreArray(da, g, &g_vec);
}

void Zeropad_Function(BGY3dH2OData BHD, Vec g, real ZP, real shift)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], dim, border, N[3];
  PetscScalar ***g_vec;
  real h[3], interval[2]; //r[3], r_s, s, r_rl_2, rr_rl_2r, rr_rl_3;

  PD = BHD->PD;
  da = BHD->da;

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];

  interval[0] = PD->interval[0];

  border = (int) ceil( ((PD->interval[1]-PD->interval[0])-(2.*ZP))/PD->h[0]/2. );

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, g, &g_vec);
   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  if( ( i[0]<=border+1 || i[0]>=N[0]-1-border) ||
	      ( i[1]<=border+1 || i[1]>=N[1]-1-border) ||
	      ( i[2]<=border+1 || i[2]>=N[2]-1-border) )
	    g_vec[i[2]][i[1]][i[0]] = shift;

/* 	  FOR_DIM */
/* 	    r[dim] = i[dim]*h[dim]+interval[0]; */
/* 	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) ); */

/* 	  if( r[0] < ZP && r[0] >= -ZP && */
/* 	      r[1] < ZP && r[1] >= -ZP && */
/* 	      r[2] < ZP && r[2] >= -ZP ) */
/* 	    g_vec[i[2]][i[1]][i[0]] = g_vec[i[2]][i[1]][i[0]]; */
/* 	  else  */
/* 	    g_vec[i[2]][i[1]][i[0]] = shift; */
	}
  DAVecRestoreArray(da, g, &g_vec);
}

// void ComputeFFTfromCoulomb(BGY3dH2OData BHD, Vec uc, Vec f_l[3],
// 			   fftw_complex *fft_data,
// 			   void *LJ_params, real damp)
void ComputeFFTfromCoulomb(BGY3dH2OData BHD, Vec uc, Vec f_l[3],
			   fftw_complex *fft_data,
			   real q2, real damp)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], ic[3], N[3], dim, index;
  PetscScalar ***v_vec;
  // real r[3], r_s, h[3], interval[2], k, fac, L, q2, sign, h3;
  real r[3], r_s, h[3], interval[2], k, fac, L, sign, h3;
  fftw_complex *tmp_fft, *(fg_fft[3]);

  PD = BHD->PD;
  da = BHD->da;
  tmp_fft = BHD->g_fft;
  FOR_DIM
    fg_fft[dim] = BHD->fg2_fft[dim];
  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  h3 = h[0]*h[1]*h[2];

  // q2 = ((double*)LJ_params)[2];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  DAVecGetArray(da, uc, &v_vec);
  index=0;
   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0];


	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

/* 	  if( r[0] < ZEROPAD && r[0] >= -ZEROPAD && */
/* 	      r[1] < ZEROPAD && r[1] >= -ZEROPAD && */
/* 	      r[2] < ZEROPAD && r[2] >= -ZEROPAD ) */
	    v_vec[i[2]][i[1]][i[0]] =
	      // damp * Coulomb_long( r_s, LJ_params);
	      damp * Coulomb_long( r_s, q2);
/* 	  else */
/* 	    v_vec[i[2]][i[1]][i[0]] = 0.0; */

	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];

	    }
	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      fft_data[index].re = 0;
	      fft_data[index].im = 0;
	      FOR_DIM
		{
		  fg_fft[dim][index].re =0;
		  fg_fft[dim][index].im = 0;
		}
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]))/SQR(L);
	      fac = EPSILON0INV/M_PI/k;
	      sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);

	      /* potential */
	      fft_data[index].re =   damp * q2 * sign *fac * exp(-k*SQR(M_PI)/SQR(G));
	      fft_data[index].im = 0;


	      /* force */
	      FOR_DIM
		{
		  fg_fft[dim][index].re =0;
		  fg_fft[dim][index].im = 2.*M_PI*ic[dim]/L*fft_data[index].re;
		}


	    }
	  index++;

	}
  DAVecRestoreArray(da, uc, &v_vec);

/*   VecSum(uc, &fac); */
/*   VecShift(uc, -fac/N[0]/N[1]/N[2]); */
  //Zeropad_Function(BHD, uc, ZEROPAD, 0.0);

  //VecSum(uc, &fac);
  //PetscPrintf(PETSC_COMM_WORLD, "int(uc)= %e\n", fac);
/*   ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, uc, tmp_fft,  */
/* 			 BHD->fft_scratch, x, n, 0); */
/*   tmp_fft[0].re=0; */

  /* Copy fft_data to temporary array for back trafo */
  for(i[0]=0; i[0]<n[0]*n[1]*n[2]; i[0]++)
    {
      tmp_fft[i[0]].re=fft_data[i[0]].re;
      tmp_fft[i[0]].im=fft_data[i[0]].im;
      //fft_data[i[0]].re= h3 * tmp_fft[i[0]].re;
      //fft_data[i[0]].im= h3 * tmp_fft[i[0]].im;

      //fprintf(stderr,"%e\n", tmp_fft[i[0]].re);
    }

  /* FFT potential */
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, uc, tmp_fft,
			 BHD->fft_scratch, x, n, 0);
  VecScale(uc, 1./L/L/L);
  //VecScale(uc, 1./N[0]/N[1]/N[2]);



  /* FFT force */
  if( f_l != PETSC_NULL)
    {
      FOR_DIM
	{
	  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, f_l[dim], fg_fft[dim],
				 BHD->fft_scratch, x, n, 0);
	  VecScale(f_l[dim], 1./L/L/L);
	}
    }

/*  VecView(uc,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

//#define SPHERE_G 1.8
//#define SPHERE_R 1.6
//#define C_G 1.8

#define SPHERE_G 2.0
#define SPHERE_R 1.0
#define C_G 1.8
// void ComputeFFTfromCoulombII(BGY3dH2OData BHD, Vec f[3] , Vec f_l[3],
// 			     fftw_complex *fft_data,
// 			     void *LJ_params, real damp)
void ComputeFFTfromCoulombII(BGY3dH2OData BHD, Vec f[3] , Vec f_l[3],
			     fftw_complex *fft_data,
			     real q2, real damp)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], ic[3], N[3], dim, index;
  // real h[3], interval[2], k, fac, L, q2, sign, h3;
  real h[3], interval[2], k, fac, L, sign, h3;
  fftw_complex *(fs_fft[3]),*(fl_fft[3]);
  real fac1, fac2, fac3, fft_s, fft_l;


  PD = BHD->PD;
  da = BHD->da;

  FOR_DIM
    fs_fft[dim] = BHD->fO_fft[dim];
  FOR_DIM
    fl_fft[dim] = BHD->fH_fft[dim];
  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  h3 = h[0]*h[1]*h[2];

  // q2 = ((double*)LJ_params)[2];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));




  fac1 = 2./((2./SQR(SPHERE_G))+(4.*SQR(SPHERE_R)));

  index=0;
   /* loop over local portion of grid */
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
	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      fft_data[index].re = 0;
	      fft_data[index].im = 0;
	      FOR_DIM
		{
		  fl_fft[dim][index].re =0;
		  fl_fft[dim][index].im =0;
		  fs_fft[dim][index].re =0;
		  fs_fft[dim][index].im =0;
		}
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]))/SQR(L);
	      fac = EPSILON0INV/k/M_PI;
	      k = sqrt(k);
	      sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);
	      fac2 = cos(2.*M_PI*SPHERE_R*k)/SQR(SPHERE_G);
	      fac3 = SPHERE_R*sin(2.*M_PI*SPHERE_R*k)/M_PI/k;

	      fft_l = damp * q2 * fac * sign *
		fac1 * exp(-SQR(M_PI*k)/SQR(SPHERE_G)) * (fac2 + fac3);
	      fft_s = damp * q2 * fac * sign * exp(-SQR(k*M_PI)/SQR(C_G)) - fft_l;

/* 	      fft_data[index].re =   damp * q2 * fac * sign *  */
/* 		( 0*exp(-SQR(k*M_PI)/SQR(C_G)) + */
/* 		   fac1 * exp(-SQR(M_PI*k)/SQR(SPHERE_G)) * (fac2 + fac3)); */
	      fft_data[index].re = fft_l;
	      fft_data[index].im = 0;



	      /* force */
	      FOR_DIM
		{
		  fl_fft[dim][index].re =0;
		  fl_fft[dim][index].im = 2.*M_PI*ic[dim]/L*fft_l;

		  fs_fft[dim][index].re =0;
		  fs_fft[dim][index].im = 2.*M_PI*ic[dim]/L*fft_s;

		}


	    }
	  //fprintf(stderr,"%e\n", fft_data[index]);
	  index++;

	}

  /* FFT force */
  if( f_l != PETSC_NULL && f != PETSC_NULL)
    {
      FOR_DIM
	{
	  /* f = f_s (+f_l) */
	  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, f[dim], fs_fft[dim],
				 BHD->fft_scratch, x, n, 0);
	  VecScale(f[dim], 1./L/L/L);
	  /*  f_l */
	  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, f_l[dim], fl_fft[dim],
				 BHD->fft_scratch, x, n, 0);
	  VecScale(f_l[dim], 1./L/L/L);

	  /* f = f_s + f_l */
	  VecAXPY(f[dim], 1.0, f_l[dim]);
	  //VecCopy(f_l[dim], f[dim]);
	}
    }

/*   VecView(f[0],PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

// void ComputeFFTSoluteII(BGY3dH2OData BHD, Vec ucl , Vec ucs, void *LJ_params,
// 			real damp, real zpad)
void ComputeFFTSoluteII(BGY3dH2OData BHD, Vec ucl , Vec ucs, real q2,
			real damp, real zpad)
{
  DA da;
  PData PD;
  int x[3], n[3], i[3], ic[3], N[3], dim, index;
  // real h[3], interval[2], k, fac, L, q2, sign, h3;
  real h[3], interval[2], k, fac, L, sign, h3;
  fftw_complex  *fs_fft,*fl_fft;
  real fac1, fac2, fac3, fft_s, fft_l;


  PD = BHD->PD;
  da = BHD->da;
  fs_fft = BHD->g_fft;
  fl_fft = BHD->gfg2_fft;

  FOR_DIM
    h[dim] = PD->h[dim];
  FOR_DIM
    N[dim] = PD->N[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  h3 = h[0]*h[1]*h[2];

  // q2 = ((double*)LJ_params)[2];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));




  fac1 = 2./((2./SQR(SPHERE_G))+(4.*SQR(SPHERE_R)));

  index=0;
   /* loop over local portion of grid */
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
	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      fl_fft[index].re = 0;
	      fl_fft[index].im = 0;
	      fs_fft[index].re = 0;
	      fs_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]))/SQR(L);
	      fac = EPSILON0INV/k/M_PI;
	      k = sqrt(k);
	      sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);
	      fac2 = cos(2.*M_PI*SPHERE_R*k)/SQR(SPHERE_G);
	      fac3 = SPHERE_R*sin(2.*M_PI*SPHERE_R*k)/M_PI/k;

	      fft_l = damp * q2 * fac * sign *
		fac1 * exp(-SQR(M_PI*k)/SQR(SPHERE_G)) * (fac2 + fac3);
	      fft_s = damp * q2 * fac * sign * exp(-SQR(k*M_PI)/SQR(C_G)) - fft_l;


	      fl_fft[index].re = fft_l;
	      fl_fft[index].im = 0;

	      fs_fft[index].re = fft_s;
	      fs_fft[index].im = 0;



	    }
	  //fprintf(stderr,"%e\n", fft_data[index]);
	  index++;

	}


  /* ucl */
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, ucl, fl_fft,
			 BHD->fft_scratch, x, n, 0);
  VecScale(ucl, 1./L/L/L);
  /*  ucs */
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, ucs, fs_fft,
			 BHD->fft_scratch, x, n, 0);
  VecScale(ucs, 1./L/L/L);

  VecSet(BHD->v[1], 0.0);
  ImposeLaplaceBoundary(BHD, ucl, BHD->v[0], BHD->v[1], zpad, NULL);
  VecSet(BHD->v[1], 0.0);
  ImposeLaplaceBoundary(BHD, ucs, BHD->v[0], BHD->v[1], zpad, NULL);

/*   VecView(ucs,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

void ComputeInitialGuess(BGY3dH2OData BHD, Vec dgO, Vec dgH, Vec dgHO, real damp)
{
  DA da;
  PData PD;
  PetscScalar ***dgH_vec, ***dgHO_vec, ***dgO_vec;
  real r[3], r_s, h[3], interval[2], beta;
  int x[3], n[3], i[3], dim;
  real q2H, q2O, q2HO;

  q2H = BHD->LJ_paramsH[2];
  q2O = BHD->LJ_paramsO[2];
  q2HO = BHD->LJ_paramsHO[2];

  PD = BHD->PD;
  da = BHD->da;


  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  DAVecGetArray(da, dgH, &dgH_vec);
  DAVecGetArray(da, dgO, &dgO_vec);
  DAVecGetArray(da, dgHO, &dgHO_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0];


	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );



	  /* Coulomb long */
	  dgH_vec[i[2]][i[1]][i[0]] =
	  //  -damp * beta* Coulomb_long( r_s, BHD->LJ_paramsH);
	    -damp * beta* Coulomb_long( r_s, q2H);
	  dgO_vec[i[2]][i[1]][i[0]] =
	  //  -damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO);
	    -damp * beta* Coulomb_long( r_s, q2O);
	  dgHO_vec[i[2]][i[1]][i[0]] =
	  //  -damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO);
	    -damp * beta* Coulomb_long( r_s, q2HO);

	}

  DAVecRestoreArray(da, dgH, &dgH_vec);
  DAVecRestoreArray(da, dgO, &dgO_vec);
  DAVecRestoreArray(da, dgHO, &dgHO_vec);


  VecSum(dgH, &BHD->ucH_0);
  VecSum(dgO, &BHD->ucO_0);
  VecSum(dgHO, &BHD->ucHO_0);
  BHD->ucH_0 *= -h[0]*h[1]*h[2]/pow(PD->interval[1]-PD->interval[0],3);
  BHD->ucO_0 *= -h[0]*h[1]*h[2]/pow(PD->interval[1]-PD->interval[0],3);
  BHD->ucHO_0 *= -h[0]*h[1]*h[2]/pow(PD->interval[1]-PD->interval[0],3);

  VecSet(dgH,0.0);
  VecSet(dgO,0.0);
  VecSet(dgHO,0.0);

/*   VecShift(BHD->ucH, BHD->ucH_0); */
/*   VecShift(BHD->ucO, BHD->ucO_0); */
/*   VecShift(BHD->ucHO, BHD->ucHO_0); */

/*   VecCopy(BHD->ucH, dgH); */
/*   VecScale(dgH, -damp*beta); */
/*   VecCopy(BHD->ucO, dgO); */
/*   VecScale(dgO, -damp*beta); */
/*   VecView(BHD->ucHO,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */
}


void RecomputeInitialData(BGY3dH2OData BHD, real damp, real damp_LJ)
{
  DA da;
  PData PD;
  PetscScalar ***gHini_vec, ***gOini_vec, ***gHOini_vec;
  PetscScalar ***(fH_vec[3]),***(fO_vec[3]),***(fHO_vec[3]);
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]),***(fHOl_vec[3]);
  PetscScalar ***cH_vec, ***cHO_vec, ***cO_vec;
  real r[3], r_s, h[3], interval[2], beta, L;
  int x[3], n[3], i[3], dim, k;
  // local LJ params and charge product
  real q2H, q2O, q2HO;
  real epsilonH, epsilonO, epsilonHO;
  real sigmaH, sigmaO, sigmaHO;


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

  PetscPrintf(PETSC_COMM_WORLD,"Recomputing initial data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);


  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;

  real periodic[27][3]={{0,0,0},{L,0,0},{-L,0,0},{0,L,0},{0,-L,0},{0,0,L},{0,0,-L},
		       {L,L,0},{-L,L,0},{L,-L,0},{-L,-L,0},
                       {L,L,L},{-L,L,L},{L,-L,L},{-L,-L,L},
                       {L,L,-L},{-L,L,-L},{L,-L,-L},{-L,-L,-L},
                       {0,L,L},{0,-L,L},{0,L,-L},{0,-L,-L},
                       {L,0,L},{-L,0,L},{L,0,-L},{-L,0,-L}};



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
  VecSet(BHD->gH_ini, 0.0);
  VecSet(BHD->gO_ini, 0.0);
  VecSet(BHD->gHO_ini, 0.0);

  /*********************************************/
  /* Compute fft from Coulomb potential (long) */
  /********************************************/
  //  ComputeFFTfromCoulomb(BHD, BHD->ucHO, BHD->fHO_l, BHD->ucHO_fft,
  //			BHD->LJ_paramsHO, damp);
  //  ComputeFFTfromCoulomb(BHD, BHD->ucH, BHD->fH_l, BHD->ucH_fft,
  //			BHD->LJ_paramsH, damp);
  //  ComputeFFTfromCoulomb(BHD, BHD->ucO, BHD->fO_l, BHD->ucO_fft,
  //			BHD->LJ_paramsO, damp);
  ComputeFFTfromCoulomb(BHD, BHD->ucHO, BHD->fHO_l, BHD->ucHO_fft,
			q2HO, damp);
  ComputeFFTfromCoulomb(BHD, BHD->ucH, BHD->fH_l, BHD->ucH_fft,
			q2H, damp);
  ComputeFFTfromCoulomb(BHD, BHD->ucO, BHD->fO_l, BHD->ucO_fft,
			q2O, damp);

  FOR_DIM
    {
      VecAXPY(BHD->fHO[dim], 1.0, BHD->fHO_l[dim]);
      VecAXPY(BHD->fH[dim], 1.0, BHD->fH_l[dim]);
      VecAXPY(BHD->fO[dim], 1.0, BHD->fO_l[dim]);
    }

  /**********************************************/

  DAVecGetArray(da, BHD->gH_ini, &gHini_vec);
  DAVecGetArray(da, BHD->gO_ini, &gOini_vec);
  DAVecGetArray(da, BHD->gHO_ini, &gHOini_vec);
  DAVecGetArray(da, BHD->cH, &cH_vec);
  DAVecGetArray(da, BHD->cO, &cO_vec);
  DAVecGetArray(da, BHD->cHO, &cHO_vec);
  FOR_DIM
    {
      DAVecGetArray(da, BHD->fH[dim], &(fH_vec[dim]));
      DAVecGetArray(da, BHD->fO[dim], &(fO_vec[dim]));
      DAVecGetArray(da, BHD->fHO[dim], &(fHO_vec[dim]));
      DAVecGetArray(da, BHD->fH_l[dim], &(fHl_vec[dim]));
      DAVecGetArray(da, BHD->fO_l[dim], &(fOl_vec[dim]));
      DAVecGetArray(da, BHD->fHO_l[dim], &(fHOl_vec[dim]));
    }



  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* set force vectors */

	  for(k=0; k<1; k++)
	    {

	  FOR_DIM
	    r[dim] = i[dim]*h[dim]+interval[0]+periodic[k][dim];


	  r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );


	  /* Lennard-Jones */
	  gHini_vec[i[2]][i[1]][i[0]] +=
	    // damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsH);
	    damp_LJ * beta* Lennard_Jones( r_s, epsilonH, sigmaH);
	  gOini_vec[i[2]][i[1]][i[0]] +=
	    // damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsO);
	    damp_LJ * beta* Lennard_Jones( r_s, epsilonO, sigmaO);
	  gHOini_vec[i[2]][i[1]][i[0]] +=
	    // damp_LJ * beta* Lennard_Jones( r_s, BHD->LJ_paramsHO);
	    damp_LJ * beta* Lennard_Jones( r_s, epsilonHO, sigmaHO);

	  /* Coulomb short */
	  gHini_vec[i[2]][i[1]][i[0]] +=
	    // damp*beta* Coulomb_short( r_s, BHD->LJ_paramsH);
            // Following the form in BGY3dH2OData_Pair_malloc() 
            // to pass member of void pointer
	    damp*beta* Coulomb_short( r_s, q2H);
	  gOini_vec[i[2]][i[1]][i[0]] +=
	    // damp*beta* Coulomb_short( r_s, BHD->LJ_paramsO);
	    damp*beta* Coulomb_short( r_s, q2O);
	  gHOini_vec[i[2]][i[1]][i[0]] +=
	    // damp*beta* Coulomb_short( r_s, BHD->LJ_paramsHO);
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

	  /* Zero padding for force vectors */
/* 	  if( r[0] < ZEROPAD && r[0] >= -ZEROPAD && */
/* 	      r[1] < ZEROPAD && r[1] >= -ZEROPAD && */
/* 	      r[2] < ZEROPAD && r[2] >= -ZEROPAD ) */
	    {
	   FOR_DIM
	    {
	      /* Lennard-Jones */
	      fH_vec[dim][i[2]][i[1]][i[0]] +=
		// damp_LJ * Lennard_Jones_grad( r_s, r[dim], BHD->LJ_paramsH);
		damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonH, sigmaH);
	      fO_vec[dim][i[2]][i[1]][i[0]] +=
		// damp_LJ * Lennard_Jones_grad( r_s, r[dim], BHD->LJ_paramsO);
		damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonO, sigmaO);
	      fHO_vec[dim][i[2]][i[1]][i[0]] +=
		// damp_LJ * Lennard_Jones_grad( r_s, r[dim], BHD->LJ_paramsHO);
		damp_LJ * Lennard_Jones_grad( r_s, r[dim], epsilonHO, sigmaHO);

 	      /* Coulomb short */
	      fH_vec[dim][i[2]][i[1]][i[0]] +=
		// damp * Coulomb_short_grad( r_s, r[dim], BHD->LJ_paramsH);
                // Following the form in BGY3dH2OData_Pair_malloc() 
                // to pass member of void pointer
		damp * Coulomb_short_grad( r_s, r[dim], q2H);
	      fO_vec[dim][i[2]][i[1]][i[0]] +=
		// damp * Coulomb_short_grad( r_s, r[dim], BHD->LJ_paramsO);
		damp * Coulomb_short_grad( r_s, r[dim], q2O);
	      fHO_vec[dim][i[2]][i[1]][i[0]] +=
		// damp * Coulomb_short_grad( r_s, r[dim], BHD->LJ_paramsHO);
		damp * Coulomb_short_grad( r_s, r[dim], q2HO);

	      /* Coulomb long */
/*  	      fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

	      /* Coulomb long */
/*  	      fHl_vec[dim][i[2]][i[1]][i[0]] =  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fOl_vec[dim][i[2]][i[1]][i[0]] =  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHOl_vec[dim][i[2]][i[1]][i[0]] =  */
/* 		damp * Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

	      /* Coulomb */
/* 	      fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		damp * Coulomb_grad( r_s, r[dim], BHD->LJ_paramsHO); */

	      /* deterministic correction */
	      cH_vec[i[2]][i[1]][i[0]] =
		// exp(-beta* LJ_repulsive( r_s, BHD->LJ_paramsH));
		exp(-beta* LJ_repulsive( r_s, epsilonH, sigmaH));
	      cO_vec[i[2]][i[1]][i[0]] =
		// exp(-beta* LJ_repulsive( r_s, BHD->LJ_paramsO));
		exp(-beta* LJ_repulsive( r_s, epsilonO, sigmaO));
	      cHO_vec[i[2]][i[1]][i[0]] =
		// exp(-beta* LJ_repulsive( r_s, BHD->LJ_paramsHO));
		exp(-beta* LJ_repulsive( r_s, epsilonHO, sigmaHO));


	    }
	    }
	    }
	}

  DAVecRestoreArray(da, BHD->gH_ini, &gHini_vec);
  DAVecRestoreArray(da, BHD->gO_ini, &gOini_vec);
  DAVecRestoreArray(da, BHD->gHO_ini, &gHOini_vec);
  DAVecRestoreArray(da, BHD->cH, &cH_vec);
  DAVecRestoreArray(da, BHD->cO, &cO_vec);
  DAVecRestoreArray(da, BHD->cHO, &cHO_vec);
  FOR_DIM
    {
      DAVecRestoreArray(da, BHD->fH[dim], &(fH_vec[dim]));
      DAVecRestoreArray(da, BHD->fO[dim], &(fO_vec[dim]));
      DAVecRestoreArray(da, BHD->fHO[dim], &(fHO_vec[dim]));
      DAVecRestoreArray(da, BHD->fH_l[dim], &(fHl_vec[dim]));
      DAVecRestoreArray(da, BHD->fO_l[dim], &(fOl_vec[dim]));
      DAVecRestoreArray(da, BHD->fHO_l[dim], &(fHOl_vec[dim]));
    }





/*   VecAXPY(BHD->gH_ini, damp*beta , BHD->ucH); */
/*   VecAXPY(BHD->gO_ini, damp*beta , BHD->ucO); */

/*   VecView(BHD->fHO[0],PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

}

void Compute_dg_H2O_inter(BGY3dH2OData BHD,
			  Vec f1[3], Vec f1_l[3], Vec g1a, Vec g1b,
			  fftw_complex *coul1_fft, real rho1, real shift1,
			  Vec f2[3], Vec f2_l[3], Vec g2a, Vec g2b,
			  fftw_complex *coul2_fft, real rho2, real shift2,
			  Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3];
  fftw_complex *(fg2_fft[3]), *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, sign, confac;

  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BHD->fg2_fft[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  fac = L/(2.*M_PI);  /* BHD->f ist nur grad U, nicht F=-grad U  */
  confac = SQR(M_PI/L/2.);

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* rho1*F1*g1a g1b*/
  /************************************************/

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BHD->v[dim], g1a, f1[dim]);
      /* special treatment: Coulomb long */
      VecAXPY(BHD->v[dim], -1.0, f1_l[dim]);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim],
			     scratch, x, n, 0);

    }
/*      VecView(f1_l[0],PETSC_VIEWER_STDERR_WORLD);   */
/*      exit(1);   */

  /* fft(g) */
  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g1b, g_fft, scratch,
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
	      dg_fft[index].re = 0;//coul1_fft[0].re*h*(g_fft[0].re-N[0]*N[1]*N[2]);
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



	      /*****************************/
	      /* long range Coulomb part */
	      dg_fft[index].re += h*sign*
		(coul1_fft[index].re*g_fft[index].re
		 - coul1_fft[index].im*g_fft[index].im);

	      dg_fft[index].im += h*sign*
		(coul1_fft[index].re*g_fft[index].im
		 + coul1_fft[index].im*g_fft[index].re);
	      /******************************/

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5))  */
/*  		{  */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].re);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft,
			 scratch, x, n, 0.0);

  VecScale(dg_help, rho1*PD->beta/L/L/L);

/*   VecScale(dg_help, 1./L/L/L); */
/*   Compute_dg_H2O_normalization_inter( BHD, g1a, g1b, BHD->v[0], BHD->v[1]); */
/*   VecPointwiseDivide(dg_help, dg_help, BHD->v[0]); */
/*   VecScale(dg_help, rho1*PD->beta); */

  VecCopy(dg_help, dg);

/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

  /************************************************/
  /* F2*g2a g2b*/
  /************************************************/
  //VecCopy(g2a, dg_help);
  //ShiftVec(da, dg_help, BHD->v[0], N);

  /* fft(f*g) */
  FOR_DIM
    {
      VecPointwiseMult(BHD->v[dim], g2a, f2[dim]);
      /* special treatment: Coulomb long */
      VecAXPY(BHD->v[dim], -1.0, f2_l[dim]);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }

  /* fft(g-1) */
  //VecCopy(g2b, dg_help);
  //VecShift(dg_help, -1.0);
  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g2b, g_fft, scratch,
			 x, n, 0);
/*   VecView(BHD->v[0],PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

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
	      dg_fft[index].re = 0;//coul2_fft[0].re*h*(g_fft[0].re-N[0]*N[1]*N[2]);
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
		 + fg2_fft[dim][index].im * g_fft[index].re);

	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac * sign *
		(-fg2_fft[dim][index].re * g_fft[index].re
		 + fg2_fft[dim][index].im * g_fft[index].im);

	      /*****************************/
	      /* long range Coulomb part */
	      dg_fft[index].re += h*sign*
		(coul2_fft[index].re*(g_fft[index].re + 0*sign/(h*rho2))
		 - coul2_fft[index].im*g_fft[index].im);

	      dg_fft[index].im += h*sign*
		(coul2_fft[index].re*g_fft[index].im
		 + coul2_fft[index].im*(g_fft[index].re + 0*sign/(h*rho2)));
	      /******************************/
	      /*****************************/
	      /* long range Coulomb part */
/* 	      dg_fft[index].re += h*sign* */
/* 		(coul2_fft[index].re*(g_fft[index].re + sign/(h*rho2)*exp(-k*confac)) */
/* 		 - coul2_fft[index].im*g_fft[index].im); */

/* 	      dg_fft[index].im += h*sign*  */
/* 		(coul2_fft[index].re*g_fft[index].im */
/* 		 + coul2_fft[index].im*(g_fft[index].re + sign/(h*rho2)*exp(-k*confac))); */
	      /******************************/

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5))  */
/*  		{  */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  //fprintf(stderr,"%e\n",dg_fft[index].re);
	  index++;
	}


  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, rho2*PD->beta/L/L/L);

/*   VecScale(dg_help, 1./L/L/L); */
/*   Compute_dg_H2O_normalization_inter( BHD, g2a, g2b, BHD->v[0], BHD->v[1]); */
/*   VecPointwiseDivide(dg_help, dg_help, BHD->v[0]); */
/*   VecScale(dg_help, rho2*PD->beta); */

/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */
  //VecShift(dg_help, PD->beta*rho2*shift2);

  VecAXPY(dg,1.0, dg_help);


/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1);  */

}






/* Compute intramolecular part */
void Compute_dg_H2O_intra(BGY3dH2OData BHD, Vec f[3], Vec f_l[3], Vec g1, Vec g2,
			  fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3];
  fftw_complex *(fg2_fft[3]), *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, beta;




  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BHD->fg2_fft[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/
  //VecCopy(g1, dg_help);
  //ShiftVec(da, dg_help, BHD->v[0], N);

  /* fft(f1*g1) */
  FOR_DIM
    {
      VecPointwiseMult(BHD->v[dim], g1, f[dim]);
      /* special treatment: Coulomb long */
      VecAXPY(BHD->v[dim], -1.0, f_l[dim]);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }

  /* fft(g2) */

/*   ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g2, g_fft, scratch, */
/* 			 x, n, 0); */

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
	      dg_fft[index].re = 0; //-coul_fft[0].re*h;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;





 	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac *
		(h*fg2_fft[dim][index].im * sin(k)/k);

	      /* - should be correct ??? */
	      //dg_fft[index].re -= h*g_fft[index].re*sin(k)/k;



	      FOR_DIM
		dg_fft[index].im += ic[dim] * k_fac *
		(-h*fg2_fft[dim][index].re * sin(k)/k);


	      /*****************************/
	      /* long range Coulomb part */
	      dg_fft[index].re +=
		coul_fft[index].re*sin(k)/k;

	      dg_fft[index].im +=
		coul_fft[index].im*sin(k)/k;
	      /******************************/

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5))  */
/*  		{  */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

		}
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
			 x, n, 0.0);

  VecScale(dg_help, PD->beta/L/L/L);


  //VecAXPY(dg, 1.0, dg_help);
  VecCopy(dg_help,dg);

  //ShiftVec(da, dg, BHD->v[0], N);
/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}


/* Compute intramolecular part */
void Compute_dg_H2O_intraII(BGY3dH2OData BHD, Vec f[3], Vec f_l[3], Vec g1, Vec tg,
			    fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3]; //, local_size;
  fftw_complex *(fg2_fft[3]), *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, beta;
  // PetscScalar *v_vec, *tg_vec;


  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BHD->fg2_fft[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */

  //PetscPrintf(PETSC_COMM_WORLD, "ACHTUNG!!: Reiner Coulomb part ist falsch!\n");
  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/

  /* fft(f1*g1) */
  FOR_DIM
    {
      VecPointwiseMult(BHD->v[dim], g1, f[dim]);
      /* special treatment: Coulomb long */
      VecAXPY(BHD->v[dim], -1.0, f_l[dim]);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }


  /* Compute int(F*g*omega) */
  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  dg_fft[index].re=0;
	  dg_fft[index].im=0;
	  FOR_DIM
	    {
	      if( i[dim] <= N[dim]/2)
		ic[dim] = i[dim];
	      else
		ic[dim] = i[dim] - N[dim];
	    }

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      FOR_DIM
		fg2_fft[dim][index].re =0 ;
	      FOR_DIM
		fg2_fft[dim][index].im =0 ;
	      dg_fft[index].re = 0;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k = 2.0*M_PI*sqrt(k)*rab/L;

	      FOR_DIM
		fg2_fft[dim][index].re = h*fg2_fft[dim][index].re * sin(k)/k;
	      FOR_DIM
		fg2_fft[dim][index].im = h*fg2_fft[dim][index].im * sin(k)/k;


	      /*****************************/
	      /* long range Coulomb part */
	      dg_fft[index].re =
		coul_fft[index].re*sin(k)/k;

	      dg_fft[index].im =
		coul_fft[index].im*sin(k)/k;
	      /******************************/

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		  FOR_DIM */
/* 		    { */
/* 		      fg2_fft[dim][index].re=0; */
/* 		      fg2_fft[dim][index].im=0; */
/* 		    } */
/* 		} */

	    }
	  //fprintf(stderr,"%e\n", dg_fft[index].re);
	  index++;
	}

  FOR_DIM
    {
      ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0.0);
      VecScale(BHD->v[dim], 1./L/L/L);
    }

  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg, dg_fft, scratch,
			 x, n, 0.0);
  VecScale(dg, PD->beta/L/L/L);
  //VecSet(dg, 0.0);

  /* int(..)/tg */
  /* Laplace^-1 * divergence */


  /*******************************************/
  /*    int/tg und dg/tg */
  /*******************************************/

/*   VecGetLocalSize(dg, &local_size); */
/*   VecGetArray( dg , &v_vec); */
/*   VecGetArray( tg, &tg_vec); */
/*   for(index=0; index<local_size; index++) */
/*     { */
/*       k = tg_vec[index];  */
/*       if( v_vec[index]<0){ */
/* 	if( fabs(k) <= 0.1) */
/* 	  v_vec[index] = v_vec[index] / 0.1; */
/* 	else */
/* 	  v_vec[index] = v_vec[index] / k; */

/*       } */
/*       else */
/* 	{ */
/* 	  if( fabs(k) <= 1.0e-4) */
/* 	    v_vec[index] = v_vec[index] / 1.0e-4 ; */
/* 	  else */
/* 	    v_vec[index] = v_vec[index] / k; */
/* 	  if( v_vec[index] > 50.) */
/* 	    v_vec[index]=50.0; */
/* 	} */

/*     } */
/*   VecRestoreArray( dg, &v_vec); */

  //VecPointwiseDivide( dg, dg, tg);
  FOR_DIM
    {
 /*      VecGetLocalSize(BHD->v[dim], &local_size); */
/*       VecGetArray( BHD->v[dim] , &v_vec); */


/*       for(index=0; index<local_size; index++) */
/* 	{ */
/*  	  k = tg_vec[index];  */
/* 	  if( v_vec[index]<0) */
/* 	    { */
/* 	      if( fabs(k) <= 0.1) */
/* 		v_vec[index] = v_vec[index] / 0.1 ; */
/* 	      else */
/* 		v_vec[index] = v_vec[index] / k; */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      if( fabs(k) <= 1.0e-4) */
/* 		v_vec[index] = v_vec[index] / 1.0e-4 ; */
/* 	      else */
/* 		v_vec[index] = v_vec[index] / k; */
/* 	      if( v_vec[index] > 50.) */
/* 		v_vec[index]=50.0; */
/* 	    } */

/* 	} */
/*       VecRestoreArray( BHD->v[dim], &v_vec); */


      VecPointwiseDivide( BHD->v[dim], BHD->v[dim], tg);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }
/*   VecRestoreArray( tg, &tg_vec); */

  /*******************************************/

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
	      k_fac = fac/k;


 	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac * h * fg2_fft[dim][index].im;

	      FOR_DIM
		dg_fft[index].im -= ic[dim] * k_fac * h * fg2_fft[dim][index].re;

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
			 x, n, 0.0);

  VecScale(dg_help, PD->beta/L/L/L);


  //VecCopy(dg_help,dg);
  VecAXPY(dg, 1.0, dg_help);

/*   VecView(dg_help ,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

/* Compute intramolecular part */
void Compute_dg_H2O_intraIII(BGY3dH2OData BHD, Vec f[3], Vec f_l[3], Vec g1, Vec tg,
			    fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3], local_size;
  fftw_complex *(fg2_fft[3]), *g_fft, *dg_fft, *scratch;
  real fac, k_fac, L, k, h, beta;
  PetscScalar *v_vec, *tg_vec;


  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];
  FOR_DIM
    fg2_fft[dim] = BHD->fg2_fft[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/

  /* fft(f1*g1) */
  FOR_DIM
    {
      VecPointwiseMult(BHD->v[dim], g1, f[dim]);
      /* special treatment: Coulomb long */
      VecAXPY(BHD->v[dim], -1.0, f_l[dim]);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }


  /* Compute int(F*g*omega) */
  index=0;
  /* loop over local portion of grid */
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

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      FOR_DIM
		fg2_fft[dim][index].re =0 ;
	      FOR_DIM
		fg2_fft[dim][index].im =0 ;

	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k = 2.0*M_PI*sqrt(k)*rab/L;

	      FOR_DIM
		fg2_fft[dim][index].re = h*fg2_fft[dim][index].re * sin(k)/k;
	      FOR_DIM
		fg2_fft[dim][index].im = h*fg2_fft[dim][index].im * sin(k)/k;




/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */

/* 		  FOR_DIM */
/* 		    { */
/* 		      fg2_fft[dim][index].re=0; */
/* 		      fg2_fft[dim][index].im=0; */
/* 		    } */
/* 		} */

	    }
	  //fprintf(stderr,"%e\n", dg_fft[index].re);
	  index++;
	}

  FOR_DIM
    {
      ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0.0);
      VecScale(BHD->v[dim], 1./L/L/L);
    }




  /* int(..)/tg */
  /* Laplace^-1 * divergence */


  /*******************************************/
  /*    int/tg  */
  /*******************************************/
  VecGetArray( tg, &tg_vec);
  FOR_DIM
    {
      VecGetLocalSize(BHD->v[dim], &local_size);
      VecGetArray( BHD->v[dim] , &v_vec);


      for(index=0; index<local_size; index++)
	{
 	  k = tg_vec[index];
	  if( k<NORM_REG)
	    v_vec[index] = v_vec[index] / NORM_REG;
	  else v_vec[index] = v_vec[index] / k;

/* 	  if( v_vec[index]<0) */
/* 	    { */
/* 	      if( fabs(k) <= 0.8) */
/* 		v_vec[index] = v_vec[index] / 0.8; */
/* 	      else */
/* 		v_vec[index] = v_vec[index] / k;	 */

/* 	    } */

/* 	  else */
/* 	    { */
/* 	      if( fabs(k) <= 1.0e-4) */
/* 		v_vec[index] = v_vec[index] / 1.0e-4 ; */
/* 	      else */
/* 		v_vec[index] = v_vec[index] / k; */
/* 	      if( v_vec[index] > 20.) */
/* 		v_vec[index]=20.0; */
/* 	    } */

	}
      VecRestoreArray( BHD->v[dim], &v_vec);

      //VecPointwiseDivide( BHD->v[dim], BHD->v[dim], tg);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }
   VecRestoreArray( tg, &tg_vec);

  /*******************************************/

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
	      FOR_DIM
		{
		  fg2_fft[dim][index].im=0;
		  fg2_fft[dim][index].re=0;
		}
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;

 	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac * h * fg2_fft[dim][index].im;

	      FOR_DIM
		dg_fft[index].im -= ic[dim] * k_fac * h * fg2_fft[dim][index].re;


	      FOR_DIM
		fg2_fft[dim][index].re = -ic[dim]/fac * coul_fft[index].im*sin(k)/k;
	      FOR_DIM
		fg2_fft[dim][index].im =  ic[dim]/fac * coul_fft[index].re*sin(k)/k;;


/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		  FOR_DIM */
/* 		    { */
/* 		      fg2_fft[dim][index].re=0; */
/* 		      fg2_fft[dim][index].im=0; */
/* 		    } */
/* 		} */

	    }
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
			 x, n, 0.0);

  VecScale(dg_help, PD->beta/L/L/L);

  /* back transformation of coulomb part, divide by tg and forward transfromation */
  VecGetArray( tg, &tg_vec);
  FOR_DIM
    {
      ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0.0);
      VecScale(BHD->v[dim], 1./L/L/L);

     VecGetLocalSize(BHD->v[dim], &local_size);
     VecGetArray( BHD->v[dim] , &v_vec);


      for(index=0; index<local_size; index++)
	{

	  k = tg_vec[index];
	  if( k<NORM_REG)
	    v_vec[index] = v_vec[index] / NORM_REG;
	  else v_vec[index] = v_vec[index] / k;
/* 	  if( v_vec[index]<0) */
/* 	    { */
/* 	      if( fabs(k) <= 0.8) */
/* 		v_vec[index] = v_vec[index] / 0.8 ; */
/* 	      else */
/* 		v_vec[index] = v_vec[index] / k; */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      if( fabs(k) <= 1.0e-4) */
/* 		v_vec[index] = v_vec[index] / 1.0e-4 ; */
/* 	      else */
/* 		v_vec[index] = v_vec[index] / k; */
/* 	      if( v_vec[index] > 20.) */
/* 		v_vec[index]=20.0; */
/* 	    } */

	}
      VecRestoreArray( BHD->v[dim], &v_vec);

      //VecPointwiseDivide(BHD->v[dim], BHD->v[dim], tg);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], fg2_fft[dim], scratch,
			     x, n, 0);
    }
  VecRestoreArray( tg, &tg_vec);

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
	      k_fac = fac/k;


 	      FOR_DIM
		dg_fft[index].re += ic[dim] * k_fac * h * fg2_fft[dim][index].im;

	      FOR_DIM
		dg_fft[index].im -= ic[dim] * k_fac * h * fg2_fft[dim][index].re;




/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg, dg_fft, scratch,
			 x, n, 0.0);

  VecScale(dg, PD->beta/L/L/L);


  VecAXPY(dg, 1.0, dg_help);

/*   VecView(dg_help ,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

/* Compute intramolecular part */
void Compute_dg_H2O_intra_ln(BGY3dH2OData BHD, Vec g, real rab, Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3], local_size;
  fftw_complex *g_fft, *dg_fft, *scratch;
  real fac, L, k, h, beta;
  PetscScalar *g_vec;


  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  fac = L/(2.*M_PI); /* siehe oben ... */
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* F(g) */
  /************************************************/


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
	      dg_fft[index].re = g_fft[0].re*h;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      //k_fac = beta*fac/k;
	      k = 2.0*M_PI*rab/L*sqrt(k);

	      /* + should be correct ??? */
	      dg_fft[index].re += h*g_fft[index].re*sin(k)/k;

	      dg_fft[index].im += h*g_fft[index].im*sin(k)/k;

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  //fprintf(stderr,"%e\n",dg_fft[index].re);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, 1./L/L/L);

/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */


  /* ln(g) */
  VecGetArray( dg_help, &g_vec);
  VecGetLocalSize(dg_help, &local_size);

  for(index=0; index<local_size; index++)
    {

      k= g_vec[index];
      if( k<1.0e-8)
	k=1.0e-8;

      g_vec[index] = -log(k);
      //g_vec[index] = -k;
    }
  VecRestoreArray(dg_help, &g_vec);
  /******************************/

  /* insure normalization condition int(u)=0 */
/*   VecSum(dg_help, &k); */
/*   VecShift( dg_help, -k/N[0]/N[1]/N[2]); */


  //VecAXPY(dg, 1.0, dg_help);
  VecCopy(dg_help,dg);


/*   VecView(dg_help,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

/* Compute intramolecular part */
void Compute_dg_H2O_intra_lnII(BGY3dH2OData BHD, Vec g, Vec t, real rab, Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3], local_size;
  fftw_complex *g_fft, *t_fft, *dg_fft, *scratch, *(f_fft[3]);
  real fac, L, k, h, beta, k_fac;
  PetscScalar *g_vec; //, *v_vec;


  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  t_fft = BHD->wHO_fft;
  dg_fft = BHD->gfg2_fft;
  FOR_DIM
    f_fft[dim] = BHD->fg2_fft[dim];

  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  fac = L/(2.*M_PI);
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* F(g) */
  /************************************************/


  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g, g_fft, scratch,
			 x, n, 0);
  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, t, t_fft, scratch,
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
	      dg_fft[index].re = g_fft[0].re*h;
	      dg_fft[index].im = 0;
	      FOR_DIM
		{
		  f_fft[dim][index].re=0;
		  f_fft[dim][index].im=0;
		}

	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      //k_fac = beta*fac/k;
	      k = 2.0*M_PI*rab/L*sqrt(k);


	      dg_fft[index].re = h*g_fft[index].re*sin(k)/k;
	      dg_fft[index].im = h*g_fft[index].im*sin(k)/k;

	      FOR_DIM
		f_fft[dim][index].re = - ic[dim]/fac *h*t_fft[index].im*sin(k)/k;

	      FOR_DIM
		f_fft[dim][index].im = + ic[dim]/fac *h*t_fft[index].re*sin(k)/k;


/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  //fprintf(stderr,"%e\n",dg_fft[index].re);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, 1./L/L/L);

  /* g>=0 ?? */
  VecGetArray( dg_help, &g_vec);
  VecGetLocalSize(dg_help, &local_size);

  for(index=0; index<local_size; index++)
    {

      k= g_vec[index];
      if( k<NORM_REG)
	k=NORM_REG;

      g_vec[index] = k;

    }
  VecRestoreArray(dg_help, &g_vec);
  /******************************/

  FOR_DIM
    {
      ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, BHD->v[dim], f_fft[dim], scratch,
		    x, n, 0.0);
      VecScale(BHD->v[dim], 1./L/L/L);



      VecPointwiseDivide(BHD->v[dim], BHD->v[dim], dg_help);

      ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, BHD->v[dim], f_fft[dim], scratch,
			     x, n, 0);
    }


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
	      k_fac = fac/k;


 	      FOR_DIM
		dg_fft[index].re +=  ic[dim] * k_fac * h * f_fft[dim][index].im;

	      FOR_DIM
		dg_fft[index].im -=  ic[dim] * k_fac * h * f_fft[dim][index].re;




/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  index++;
	}

  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg, dg_fft, scratch,
			 x, n, 0.0);

  /* Vorzeichen direkt aus BGY3dM equation */
  VecScale(dg, -1./L/L/L);



/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}

/* Compute intramolecular part */
void Compute_dg_H2O_intra_lnIII(BGY3dH2OData BHD, Vec g, Vec t, real rab, Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3], local_size;
  fftw_complex *g_fft, *t_fft, *dg_fft, *scratch, *(f_fft[3]), *n_fft;
  real fac, L, k, h, beta, k_fac;
  PetscScalar *g_vec; //, *v_vec;


  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  t_fft = BHD->wHO_fft;
  dg_fft = BHD->gfg2_fft;
  n_fft = BHD->fg2_fft[0];

  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  fac = L/(2.*M_PI);
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /************************************************/
  /* F(g) */
  /************************************************/


  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g, g_fft, scratch,
			 x, n, 0);
  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, t, t_fft, scratch,
			 x, n, 0);

  index=0;
  /* loop over local portion of grid */
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

	  if( ic[0]==0 && ic[1]==0 && ic[2]==0)
	    {
	      dg_fft[index].re = t_fft[0].re*h;
	      dg_fft[index].im = 0;

	      n_fft[index].re = g_fft[0].re*h;
	      n_fft[index].im = 0;


	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      //k_fac = beta*fac/k;
	      k = 2.0*M_PI*rab/L*sqrt(k);


	      n_fft[index].re = h*g_fft[index].re*sin(k)/k;
	      n_fft[index].im = h*g_fft[index].im*sin(k)/k;

	      g_fft[index].re = h*t_fft[index].re*sin(k)/k;
	      g_fft[index].im = h*t_fft[index].im*sin(k)/k;




	    }
	  //fprintf(stderr,"%e\n",dg_fft[index].re);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, 1./L/L/L);

  /* g>=0 ?? */
  VecGetArray( dg_help, &g_vec);
  VecGetLocalSize(dg_help, &local_size);

  for(index=0; index<local_size; index++)
    {

      k= g_vec[index];
      if( k<NORM_REG)
	k=NORM_REG;

      g_vec[index] = k;

    }
  VecRestoreArray(dg_help, &g_vec);
  /******************************/

  FOR_DIM
    {
      ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, BHD->v[dim], f_fft[dim], scratch,
		    x, n, 0.0);
      VecScale(BHD->v[dim], 1./L/L/L);

    }


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
	      FOR_DIM
		{
		  f_fft[dim][index].re=0;
		  f_fft[dim][index].im=0;
		}

	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      k_fac = k*4.*SQR(M_PI)/SQR(L);
	      k = 2.0*M_PI*rab/L*sqrt(k);

	      dg_fft[index].re = -k_fac*h*t_fft[index].re*sin(k)/k;
	      dg_fft[index].im = -k_fac*h*t_fft[index].im*sin(k)/k;


 	      FOR_DIM
		f_fft[dim][index].re = - ic[dim]/fac *h*g_fft[index].im*sin(k)/k;

	      FOR_DIM
		f_fft[dim][index].im = + ic[dim]/fac *h*g_fft[index].re*sin(k)/k;;




/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5)) */
/* 		{ */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */

	    }
	  index++;
	}
  FOR_DIM
    {
      ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg, f_fft[dim], scratch,
			     x, n, 0.0);
      VecScale(dg, 1./L/L/L);

      VecPointwiseMult(BHD->v[dim], BHD->v[dim], dg);

    }
  VecCopy(BHD->v[0], dg);
  VecAXPY(dg, 1.0, BHD->v[1]);
  VecAXPY(dg, 1.0, BHD->v[2]);
  VecPointwiseDivide(dg, dg, dg_help);

  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, BHD->v[0], dg_fft, scratch,
			 x, n, 0.0);
  VecScale(BHD->v[0], 1./L/L/L);

  VecAXPY(dg, 1.0, BHD->v[0]);

  VecPointwiseDivide(dg, dg, dg_help);

/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1); */

  /* inverse Laplace */
  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, dg, g_fft, scratch,
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
	      k_fac = k*4.*SQR(M_PI)/SQR(L);


	      dg_fft[index].re = -h/k_fac*g_fft[index].re;
	      dg_fft[index].im = -h/k_fac*g_fft[index].im;


	    }
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg, dg_fft, scratch,
			 x, n, 0.0);
  /* Vorzeichen direkt aus BGY3dM equation */
  VecScale(dg, -1./L/L/L);



/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

}


void Compute_Zero_Correction(BGY3dH2OData BHD, Vec dg)
{
  PData PD;
  DA da;
  int x[3], n[3], dim, N[3];
  fftw_complex *g_fft, *scratch;
  real L, h;


  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, dg, g_fft, scratch,
			 x, n, 0);
  //PetscPrintf(PETSC_COMM_WORLD,"%e\n",g_fft[0].re);
  g_fft[0].re=0;
  g_fft[0].im=0;

  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg, g_fft, scratch,
		    x, n, 0.0);
  VecScale(dg, h/L/L/L);
}

/* Compute normalization condition */
void Compute_dg_H2O_normalization_intra(BGY3dH2OData BHD, Vec g, real rab,
					Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3];
  fftw_complex  *g_fft, *dg_fft, *scratch;
  real fac, L, k, h, beta;
  PetscScalar *g_vec;
  int local_size;
  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];


  h=PD->h[0]*PD->h[1]*PD->h[2];
  g_fft = BHD->g_fft;
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* fft(g/t) */
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

	      dg_fft[index].re = h*g_fft[0].re;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      k = (SQR(ic[2])+SQR(ic[1])+SQR(ic[0]));
	      //k_fac = beta*fac/k;
	      k = 2.0*M_PI*sqrt(k)*rab/L;



	      /* + hier richtig !! */
	      dg_fft[index].re += h*g_fft[index].re*sin(k)/k;

	      dg_fft[index].im += h*g_fft[index].im*sin(k)/k;

/* 	      if( (SQR(ic[0])+SQR(ic[1])+SQR(ic[2]))>SQR(N[0]/2-5))  */
/*  		{  */
/* 		  dg_fft[index].re= 0; */
/* 		  dg_fft[index].im= 0; */
/* 		} */


	    }
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);


  VecScale(dg_help, 1./L/L/L);

  /* g>=0 ?? */
  VecGetArray( dg_help, &g_vec);
  VecGetLocalSize(dg_help, &local_size);

  for(index=0; index<local_size; index++)
    {

      k= g_vec[index];
      if( k<1.0e-8)
	k=1.0e-8;

      g_vec[index] = k;

    }
  VecRestoreArray(dg_help, &g_vec);
  /******************************/

  //VecAXPY(dg, 1.0, dg_help);
  VecCopy(dg_help,dg);


  //VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);

}

/* Compute normalization condition */
void Compute_dg_H2O_normalization_inter(BGY3dH2OData BHD, Vec g1, Vec g2,
					Vec dg, Vec dg_help)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], dim, index, N[3], ic[3];
  fftw_complex  *g1_fft, *g2_fft, *dg_fft, *scratch;
  real fac, L, h, beta, sign; // k
  // PetscScalar *g_vec;
  // int local_size;
  PD=BHD->PD;

  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];


  h=PD->h[0]*PD->h[1]*PD->h[2];
  g1_fft = BHD->g_fft;
  g2_fft = BHD->fg2_fft[0];
  dg_fft = BHD->gfg2_fft;
  scratch = BHD->fft_scratch;
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;
  fac = L/(2.*M_PI); /* siehe oben ... */


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));


  /* fft(g) */
  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g1, g1_fft, scratch,
			 x, n, 0);

  ComputeFFTfromVec_fftw(da, BHD->fft_plan_fw, g2, g2_fft, scratch,
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

	      dg_fft[index].re = 0;//h*h*g1_fft[0].re*g2_fft[0].re;
	      dg_fft[index].im = 0;
	    }
	  else
	    {
	      sign = COSSIGN(ic[0])*COSSIGN(ic[1])*COSSIGN(ic[2]);


	      /* + hier richtig !! */
	      dg_fft[index].re += h*h* sign*(
					g1_fft[index].re*g2_fft[index].re -
					g1_fft[index].im*g2_fft[index].im
					);
	      dg_fft[index].im += h*h* sign*(
					g1_fft[index].re*g2_fft[index].im +
					g1_fft[index].im*g2_fft[index].re
					);




	    }
	  //fprintf(stderr,"%e\n",fg2_fft[0][index].im);
	  index++;
	}
  ComputeVecfromFFT_fftw(da, BHD->fft_plan_bw, dg_help, dg_fft, scratch,
		    x, n, 0.0);

  VecScale(dg_help, 100./pow(L,6.0));
  VecShift(dg_help, 1.0);


  /* g>=0 ?? */
/*   VecGetArray( dg_help, &g_vec); */
/*   VecGetLocalSize(dg_help, &local_size); */

/*   for(index=0; index<local_size; index++) */
/*     { */

/*       k= g_vec[index]; */
/*       if( k<1.0e-8)  */
/* 	k=1.0e-8; */

/*       g_vec[index] = k; */

/*     } */
/*   VecRestoreArray(dg_help, &g_vec); */
  /******************************/

  //VecAXPY(dg, 1.0, dg_help);
  VecCopy(dg_help,dg);


  //VecView(dg_help,PETSC_VIEWER_STDERR_WORLD);
  //exit(1);

}

#define APP_NORM 1.0e-2
void Solve_NormalizationH2O(BGY3dH2OData BHD, Vec gH, Vec gO, Vec gHO, Vec gOH,
			 Vec tH, Vec tO, Vec tHO, Vec tOH, Vec dg, Vec dg_help)
{
  real normH, normO, normHO;
  int count=0;

  VecCopy(gH, tH);
  VecCopy(gO, tO);
  VecCopy(gHO, tHO);


  do
    {
      count++;
      Compute_dg_H2O_normalization_intra( BHD, tO, r_HO, dg, dg_help);
      VecCopy(tHO, tOH);
      VecPointwiseDivide(tHO, gHO, dg);
      VecAXPY(tOH, -1.0, tHO);
      VecNorm(tOH, NORM_INFINITY, &normHO);
      //PetscPrintf(PETSC_COMM_WORLD,"%e ", normHO);

/*       VecView(tHO,PETSC_VIEWER_STDERR_WORLD); */
/*       exit(1);  */

      Compute_dg_H2O_normalization_intra( BHD, tHO, r_HO, dg, dg_help);
      VecCopy(tO, tOH);
      VecPointwiseDivide(tO, gO, dg);
      VecAXPY(tOH, -1.0, tO);
      VecNorm(tOH, NORM_INFINITY, &normO);
      //PetscPrintf(PETSC_COMM_WORLD,"%e ", normO);

      Compute_dg_H2O_normalization_intra( BHD, tH, r_HH, dg, dg_help);
      VecCopy(tH, tOH);
      VecPointwiseDivide(tH, gH, dg);
      VecAXPY(tOH, -1.0, tH);
      VecNorm(tOH, NORM_INFINITY, &normH);
      //PetscPrintf(PETSC_COMM_WORLD,"%e \n", normH);

    }while(count<50 &&  (normH>APP_NORM || normO>APP_NORM ||normHO>APP_NORM));




  return ;
}

void Solve_NormalizationH2O_small_old(BGY3dH2OData BHD, Vec gc, real rc, Vec g, Vec t,
				  Vec dg, Vec dg_help, real zpad)
{

/*   VecCopy(g,t); */
/*   return; */

  Compute_dg_H2O_normalization_intra( BHD, gc, rc, dg, dg_help);
  //Smooth_Function(BHD, dg, zpad-1, zpad, 1.0);

  VecPointwiseDivide(t, g, dg);
  Zeropad_Function(BHD, dg, zpad, 1.0);
/*   VecView(g,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */
}

#define MAXENTRY 5.0e+1
void Solve_NormalizationH2O_small(BGY3dH2OData BHD, Vec gc, real rc, Vec g, Vec t,
				  Vec dg, Vec dg_help, real zpad)
{
  int local_size, i;
  PetscScalar *t_vec, *g_vec, *dg_vec;
  real k;

/*   VecCopy(g, t); */
/*   return; */

  Compute_dg_H2O_normalization_intra( BHD, gc, rc, dg, dg_help);
  //Smooth_Function(BHD, dg, zpad-1, zpad, 1.0);



  //VecPointwiseDivide(t, g, dg);

  VecGetLocalSize(t, &local_size);
  VecGetArray(t, &t_vec);
  VecGetArray(g, &g_vec);
  VecGetArray(dg, &dg_vec);

    for(i=0; i<local_size; i++)
    {
      k = dg_vec[i];
      if(k<NORM_REG)
	{
/* 	  if( g_vec[i]<0.1) */
/* 	    t_vec[i]=0; */
/* 	  else */
	    t_vec[i] = g_vec[i]/NORM_REG;
	}
      else
	t_vec[i] = g_vec[i]/k;

    }
    VecRestoreArray(t, &t_vec);
    VecRestoreArray(g, &g_vec);
    VecRestoreArray(dg, &dg_vec);


  //Zeropad_Function(BHD, dg, zpad, 1.0);
/*   VecView(g,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */
}


void Solve_NormalizationH2O_smallII(BGY3dH2OData BHD, Vec gc, real rc, Vec g, Vec t,
				  Vec dg, Vec dg_help, real zpad)
{
  int local_size, i;
  PetscScalar *t_vec, *g_vec, *dg_vec;
  real k;

/*   VecCopy(g, t); */
/*   return; */


  Compute_dg_H2O_normalization_intra( BHD, gc, rc, dg, dg_help);
  //Smooth_Function(BHD, dg, zpad-1, zpad, 1.0);


  //VecPointwiseDivide(t, g, dg);

  VecGetLocalSize(t, &local_size);
  VecGetArray(t, &t_vec);
  VecGetArray(g, &g_vec);
  VecGetArray(dg, &dg_vec);

  for(i=0; i<local_size; i++)
    {
      k = dg_vec[i];
      if(k<NORM_REG2)
	{
/* 	  if( g_vec[i]<0.1) */
/* 	    t_vec[i]=0; */
/* 	  else */
	  t_vec[i] = g_vec[i]/NORM_REG2;

/* 	    if( t_vec[i] >MAXENTRY) */
/* 	      t_vec[i] = MAXENTRY; */

	}
      else
	t_vec[i] = g_vec[i]/k;

    }
  VecRestoreArray(t, &t_vec);
  VecRestoreArray(g, &g_vec);
  VecRestoreArray(dg, &dg_vec);

  Zeropad_Function(BHD, dg, zpad, 1.0);
  /*   VecView(g,PETSC_VIEWER_STDERR_WORLD);  */
  /*   exit(1);  */
}

#define DAMPO 1.0
/* solve */
Vec BGY3d_solve_2site(PData PD, Vec g_ini, int vdim)
{
  BGY3dH2OData BHD;
  Vec g0H, g0O, g0HO, dgH, dgO, dgHO, dg_new, dg_new2, f, gH, gO, gHO;
  Vec dgOH, gOH;
  Vec tH, tO, tHO, tOH;
  real a=0.9, damp, damp_LJ, damp_start=0.0;
  real aH, aHO, aO, a0=0.9, a1=0.5, count=0.0, zpad; //, norm;
  int max_iter=25, iter, mycount=0, upwards=0, namecount=0; // in_iter
  char nameO[20], nameH[20], nameHO[20];
  PetscScalar dgH_norm, dgO_norm, dgHO_norm, norm_tol=1.0e-2;
  PetscScalar dgH_old, dgHO_old, dgO_old;
  real ti;
  int iteri;
  // PetscScalar dgOH_norm;
  PetscTruth kflg, load_flag;
  PetscViewer viewer;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Fourier ansatz...\n");

  PetscOptionsHasName(PETSC_NULL,"-pair",&kflg);
  if(kflg)
    {

      BHD = BGY3dH2OData_Pair_malloc(PD);
      if( r_HH >0)
	{
	  PetscPrintf(PETSC_COMM_WORLD,"WARNING: Solvent not a 2-Site model!\n");
	}
    }
  else
    {
      exit(1);
      //BHD = BGY3dFourierData_malloc(PD);
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }

  /* read BGY3d specific things from command line */
  /* Mixing parameter */
  PetscOptionsGetReal(PETSC_NULL,"-lambda",&a0, PETSC_NULL);
  if(a>1 || a<0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lambda out of range: lambda=%f\n",a);
      exit(1);
    }
   /* Get damp_start from command line*/
  PetscOptionsGetReal(PETSC_NULL,"-damp_start",&damp_start, PETSC_NULL);
  /* Number of total iterations */
  PetscOptionsGetInt(PETSC_NULL,"-max_iter",&max_iter, PETSC_NULL);
  /* norm_tol for convergence test */
  PetscOptionsGetReal(PETSC_NULL,"-norm_tol",&norm_tol, PETSC_NULL);
  /* Zeropad */
  PetscOptionsGetReal(PETSC_NULL,"-zpad",&zpad, PETSC_NULL);
  /*********************************/

  PetscPrintf(PETSC_COMM_WORLD,"lambda = %f\n",a);
  PetscPrintf(PETSC_COMM_WORLD,"tolerance = %e\n",norm_tol);
  PetscPrintf(PETSC_COMM_WORLD,"zpad = %f\n",zpad);
  PetscPrintf(PETSC_COMM_WORLD,"max_iter = %d\n",max_iter);

  DACreateGlobalVector(BHD->da, &gH);
  DACreateGlobalVector(BHD->da, &gO);
  DACreateGlobalVector(BHD->da, &gHO);
  DACreateGlobalVector(BHD->da, &dgH);
  DACreateGlobalVector(BHD->da, &dgO);
  DACreateGlobalVector(BHD->da, &dgHO);
  DACreateGlobalVector(BHD->da, &dg_new);
  DACreateGlobalVector(BHD->da, &dg_new2);
  DACreateGlobalVector(BHD->da, &f);

  DACreateGlobalVector(BHD->da, &dgOH);
  DACreateGlobalVector(BHD->da, &gOH);

  DACreateGlobalVector(BHD->da, &tH);
  DACreateGlobalVector(BHD->da, &tO);
  DACreateGlobalVector(BHD->da, &tHO);
  DACreateGlobalVector(BHD->da, &tOH);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix */
  InitializeLaplaceMatrix(BHD, zpad);
  /* Create KSP environment */
  InitializeKSPSolver(BHD);
#endif

  g0H=BHD->gH_ini;
  g0O=BHD->gO_ini;
  g0HO=BHD->gHO_ini;

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dgHO,0);
  VecSet(dgOH,0.0);


/*   VecSetRandom_H2O(dgH,1.0); */
/*   VecSetRandom_H2O(dgO,1.0); */
/*   VecSetRandom_H2O(dgHO,1.0); */


  /* load initial configuration from file ??? */
  PetscOptionsHasName(PETSC_NULL,"-loadH2O",&load_flag);
  if(load_flag)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      /* dgH */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgH.bin",
			    FILE_MODE_READ,&viewer);
      VecLoad(viewer,VECMPI, &dgH);
      PetscViewerDestroy(viewer);
      /* dgO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgO.bin",
			    FILE_MODE_READ,&viewer);
      VecLoad(viewer,VECMPI, &dgO);
      PetscViewerDestroy(viewer);
      /* dgHO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgHO.bin",
			    FILE_MODE_READ,&viewer);
      VecLoad(viewer,VECMPI, &dgHO);
      PetscViewerDestroy(viewer);
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
    }

  VecSet(dg_new,0.0);






  //ComputeH2O_g( gOH, g0HO, dgOH);


  for( damp=damp_start; damp <=1; damp+=0.1)
    {
      if(damp==-0.01)
	{
	  damp_LJ=0;
	  //a0=0.4;
	  RecomputeInitialData(BHD, 0, 1.0);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}
      else if(damp==0.0)
	{
	  damp_LJ=1.0;
	  //a0=0.5;
	  RecomputeInitialData(BHD, 0, 1.0);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}
      else
	{
	  damp_LJ=1.0;
	  //a0=0.0002/damp;
	  count+=1.0;
	  //a0 = 0.1/4./count;
	  //a0 = 0.1/count;
	  RecomputeInitialData(BHD, (damp), 1.0);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}

/*       damp=0.05; */
/*       RecomputeInitialData(BHD, damp, 1.0); */

/*       Smooth_Function(BHD, g0HO, SL, SR, 0.0); */
/*       Smooth_Function(BHD, g0O, SL, SR, 0.0); */
/*       Smooth_Function(BHD, g0H, SL, SR, 0.0); */
/*       Zeropad_Function(BHD, g0HO, zpad, 0.0); */
/*       Zeropad_Function(BHD, g0O, zpad, 0.0); */
/*       Zeropad_Function(BHD, g0H, zpad, 0.0); */
      ImposeLaplaceBoundary(BHD, g0H, tH, BHD->xH, zpad, NULL);
      ImposeLaplaceBoundary(BHD, g0O, tH, BHD->xO, zpad, NULL);
      ImposeLaplaceBoundary(BHD, g0HO, tH, BHD->xHO, zpad, NULL);
      Zeropad_Function(BHD, g0HO, zpad, 0.0);
      Zeropad_Function(BHD, g0O, zpad, 0.0);
      Zeropad_Function(BHD, g0H, zpad, 0.0);
      /* g=g0*exp(-dg) */

      ComputeH2O_g( gHO, g0HO, dgHO);
      ComputeH2O_g( gH, g0H, dgH);
      ComputeH2O_g( gO, g0O, dgO);

      a1=a0;
      a=a0;
      aH=a;
      aO=a;
      aHO=a;

  for(iter=0; iter<max_iter; iter++)
    {



      if( !(iter%10) && iter>0 )
	{
	  a=a1;
	  //PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a);
	  aH=a;
	  aO=a;
	  aHO=a;
	}
      else if( iter==20)
	{
/*  	  a=0.1;  */
/*  	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a);  */
	  aH=a;
	  aO=a;
	  aHO=a;
	}
      else
	{
	  //a=0.01;
	  //PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a);
	  a=a0;
	  aH=a;
	  aO=a;
	  aHO=a;
	}

      PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms: %e %e ", iter+1, NORM_REG, NORM_REG2);
      /* f=integral(g) */
      if(kflg)
	{
/* 	  Compute_dg_H2O_normalization_inter( BHD, gH, gH, tHO, f); */

/* 	  VecView(BHD->cHO,PETSC_VIEWER_STDERR_WORLD);          */
/*  	  exit(1);     */

	  //if(iter < 50) goto gH;
/* 	  if( !(iter%50) && iter>0 && NORM_REG>=2.0e-2) */
/* 	    NORM_REG/=2.; */

/* 	  ComputeH2O_g( gH, g0H, dgH); */
/* 	  ComputeH2O_g( gO, g0O, dgO); */
/* 	  ComputeH2O_g( gHO, g0HO, dgHO); */
	  goto gOH;
	  /* g_HO */
	  Compute_dg_H2O_inter(BHD,
			       BHD->fHO, BHD->fHO_l, gHO, gH,
			       BHD->ucHO_fft, BHD->rho_O, BHD->ucHO_0,
			       BHD->fH, BHD->fH_l, gH, gHO,
			       BHD->ucH_fft, BHD->rho_H, BHD->ucH_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);
	  //VecSet(dg_new,0.0);
	  VecPointwiseMult(dg_new, dg_new, BHD->cHO);



	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gO, tO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gO, tO, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);



#ifdef INTRA1
	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fH, BHD->fH_l, tH, PETSC_NULL,
			       BHD->ucH_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gH, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);
#endif
#ifdef INTRA2
	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gH, r_HO, tHO, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fH, BHD->fH_l, tH, tHO,
				 BHD->ucH_fft, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);
#endif

/* 	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);         */
/* 	  exit(1);    */

	  VecAXPY(dg_new, PD->beta, BHD->ucHO);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);
	  if(iter>=0)
	    {
	      ti=ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xHO, zpad, &iteri);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	      PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);
	    }

/* 	  VecNorm(dg_new, NORM_2, &norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"nrom=%e  ",norm); */

	  VecCopy(dgHO, f);
	  //VecAXPBY(dgHO, a, (1-a), dg_new);
	  VecAXPBY(dgHO, aHO, (1-aHO), dg_new);
	  VecAXPY(f, -1.0, dgHO);
	  VecNorm(f, NORM_INFINITY, &dgHO_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"HO= %e  (%f)  ",  dgHO_norm/aHO, aHO);
	  ComputeH2O_g( gHO, g0HO, dgHO);

/* 	  for(in_iter=0; in_iter<0; in_iter++) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD,"in_iter %d= ",in_iter); */

/* 	  Solve_NormalizationH2O(BHD, gH,  gO, gHO,  gOH, */
/* 				 tH, tO, tHO, tOH, dg_new2, f); */
	  goto gH;
	gOH:
	  //goto gH;
	  /* g_OH */
	  Compute_dg_H2O_inter(BHD,
			       BHD->fO, BHD->fO_l, gO, gHO,
			       BHD->ucO_fft, BHD->rho_O, BHD->ucO_0,
			       BHD->fHO, BHD->fHO_l, gHO, gH,
			       BHD->ucHO_fft, BHD->rho_H, BHD->ucHO_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);
	  VecPointwiseMult(dg_new, dg_new, BHD->cHO);


	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gH, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gH, tH, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gH, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
 	  VecAXPY(dg_new, 1.0, dg_new2);


/* 	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);         */
/* 	  exit(1);    */

#ifdef INTRA1

	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fO, BHD->fO_l, tO, PETSC_NULL,
			       BHD->ucO_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);


#endif
#ifdef INTRA2
	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gO, r_HO, tHO, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fO, BHD->fO_l, tO, tHO,
				 BHD->ucO_fft, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);
#endif

/* 	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);         */
/* 	  exit(1);   */

	  VecAXPY(dg_new, PD->beta, BHD->ucHO);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);

	  if(iter>=0)
	    {
	      ti=ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xHO, zpad, &iteri);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	      //PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);
	    }
/* 	  VecNorm(dg_new, NORM_2, &norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"nrom=%e  ",norm); */

	  VecCopy(dgHO, f);
	  //VecAXPBY(dgHO, a, (1-a), dg_new);
	  VecAXPBY(dgHO, aHO, (1-aHO), dg_new);
	  VecAXPY(f, -1.0, dgHO);
	  VecNorm(f, NORM_INFINITY, &dgHO_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"HO= %e  (%f)  ",  dgHO_norm/aHO, aHO);
	  //ComputeH2O_g( gHO, g0HO, dgHO);

	  //VecView(dgHO,PETSC_VIEWER_STDERR_WORLD);
 	  //exit(1);

	  //if(iter==max_iter-1) VecCopy(f, dgHO);
/* 	  } */
 	  /*****************************************/
/* 	  Solve_NormalizationH2O( BHD,  gH,  gO,  gHO, gOH, tH, tO, tHO, tOH, dg_new, f,  */
/* 			   norm_tol); */
	  /*****************************************/
/* 	  for(in_iter=0; in_iter<0; in_iter++) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD,"in_iter %d= ",in_iter); */
/* 	  VecCopy(dgHO, dgO); */
/* 	  VecCopy(dgHO, dgH); */
/* 	  goto ende; */

	gH:
	  //goto ende;
	  //if(damp==0.0 && iter<50) goto ende;
	  /* g_H */
  	  Compute_dg_H2O_inter(BHD,
			       BHD->fHO, BHD->fHO_l, gHO, gHO,
			       BHD->ucHO_fft, BHD->rho_O, BHD->ucHO_0,
			       BHD->fH, BHD->fH_l, gH, gH,
			       BHD->ucH_fft, BHD->rho_H, BHD->ucH_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);
	  //VecScale(dg_new, 0.5);
	  VecPointwiseMult(dg_new, dg_new, BHD->cH);


/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);        */
/* 	  exit(1);   */


	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gHO, tHO, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);




#ifdef INTRA1

 	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fHO, BHD->fHO_l, tHO, PETSC_NULL,
			       BHD->ucHO_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);

#endif
#ifdef INTRA2
	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HO, tH, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fHO, BHD->fHO_l, tHO, tH,
				 BHD->ucHO_fft, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);
#endif

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);        */
/*   	  exit(1);   */


	  VecAXPY(dg_new, PD->beta, BHD->ucH);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);


	  if(iter>=0)
	    {
	      ti=ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xH, zpad, &iteri);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	      //PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);
	    }
	  VecCopy(dgH, f);
	  //VecAXPBY(dgH, a, (1-a), dg_new);
	  VecAXPBY(dgH, aH, (1-aH), dg_new);
	  VecAXPY(f, -1.0, dgH);
	  VecNorm(f, NORM_INFINITY, &dgH_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"H= %e  (%f)  ", dgH_norm/aH, aH);
	  //ComputeH2O_g( gH, g0H, dgH);

	  //if(iter<50 )goto ende;
	  //if(iter==max_iter-1) VecCopy(f, dgH);
	  //ComputeH2O_Renormalization(BHD, gH);
/* 	  } */
/* 	  for(in_iter=0; in_iter<300; in_iter++) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD,"in_iter %d= ",in_iter); */
	gO:
	  /* g_O */
	  //goto ende;
	  Compute_dg_H2O_inter(BHD,
			       BHD->fHO, BHD->fHO_l, gHO, gHO,
			       BHD->ucHO_fft, BHD->rho_H, BHD->ucHO_0,
			       BHD->fO, BHD->fO_l, gO, gO,
			       BHD->ucO_fft, BHD->rho_O, BHD->ucO_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);
	  VecPointwiseMult(dg_new, dg_new, BHD->cO);

	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gHO, tHO, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);

/* 	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);        */
/* 	  exit(1);   */


#ifdef INTRA1

	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fHO, BHD->fHO_l, tHO, PETSC_NULL,
			       BHD->ucHO_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);

#endif
#ifdef INTRA2

	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HO, tO, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fHO, BHD->fHO_l, tHO, tO,
				 BHD->ucHO_fft, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  //Solve_NormalizationH2O_smallII( BHD, gO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);
#endif

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);       */
/* 	  exit(1);  */

	  VecAXPY(dg_new, PD->beta, BHD->ucO);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);
	  if(iter>=0)
	    {
	      ti=ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xO, zpad, &iteri);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	      //PetscPrintf(PETSC_COMM_WORLD,"%e %d ", ti, iteri);
	    }
	  VecCopy(dgO, f);
	  //VecAXPBY(dgO, a, (1-a), dg_new);
 	  VecAXPBY(dgO, aO, (1-aO), dg_new);
	  VecAXPY(f, -1.0,  dgO);
	  VecNorm(f, NORM_INFINITY, &dgO_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"O= %e  (%f)  ", dgO_norm/aO, aO);
	  //ComputeH2O_g( gO, g0O, dgO);

	ende:
	  ;
	  ComputeH2O_g( gHO, g0HO, dgHO);
	  ComputeH2O_g( gH, g0H, dgH);
	  ComputeH2O_g( gO, g0O, dgO);


/* 	  norm = ComputeCharge(BHD, gHO, gH); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " HO=%.2e ", norm); */
/* 	  norm = ComputeCharge(BHD, gO, gHO); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " O=%.2e ", norm); */
/* 	  norm = ComputeCharge(BHD, gH, gHO); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " H=%.2e ", norm); */




	}
      else {
        // nothing
      }

      //PetscPrintf(PETSC_COMM_WORLD,"\n");

      /* (fancy) step size control */
      mycount++;
      if( ((iter-1)%10) &&
	  (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
	{
	  upwards=1;
	}
      else if(iter>20 && !((iter-1)%10) && upwards==0 &&
	      (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
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
      dgH_old = dgH_norm/aH;
      dgO_old = dgO_norm/aO;
      dgHO_old = dgHO_norm/aHO;
      /*********************************/


      if(dgH_norm/aH<=norm_tol &&  dgO_norm/aO<=norm_tol && dgHO_norm/aHO<=norm_tol )//&&NORM_REG<=1.0e-2)
	break;


      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }


/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */

  /*************************************/
  /* output */
  namecount++;
  sprintf(nameH, "vecH-%d.m", namecount-1);
  sprintf(nameO, "vecO-%d.m", namecount-1);
  sprintf(nameHO, "vecHO-%d.m", namecount-1);
  PetscPrintf(PETSC_COMM_WORLD,"Writing files...");
  /* g_H */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameH,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecH.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gH,viewer);
  PetscViewerDestroy(viewer);
  /* g_b */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameO,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecO.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gO,viewer);
  PetscViewerDestroy(viewer);
  /* g_HO */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameHO,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gHO,viewer);
  PetscViewerDestroy(viewer);
  /* g_ba */
  //PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecba.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m",FILE_MODE_WRITE,&viewer);
  //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  //VecView(gba,viewer);
  //PetscViewerDestroy(viewer);
  PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  /************************************/
  /* save g to binary file */
  PetscOptionsHasName(PETSC_NULL,"-saveH2O",&load_flag);
  if(load_flag)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
      /* dgH */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgH.bin",
			    FILE_MODE_WRITE,&viewer);
      VecView(dgH,viewer);
      PetscViewerDestroy(viewer);
      /* dgO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgO.bin",
			    FILE_MODE_WRITE,&viewer);
      VecView(dgO,viewer);
      PetscViewerDestroy(viewer);
      /* dgHO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgHO.bin",
			    FILE_MODE_WRITE,&viewer);
      VecView(dgHO,viewer);
      PetscViewerDestroy(viewer);
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
    }

  /************************************/
  /* save g2 to binary file */
  PetscPrintf(PETSC_COMM_WORLD,"Writing g2 files...");
  /* g2H */
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2H.bin",
			FILE_MODE_WRITE,&viewer);
  VecView(gH,viewer);
  PetscViewerDestroy(viewer);
  /* g2O */
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2O.bin",
			FILE_MODE_WRITE,&viewer);
  VecView(gO,viewer);
  PetscViewerDestroy(viewer);
  /* g2HO */
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2HO.bin",
			FILE_MODE_WRITE,&viewer);
  VecView(gHO,viewer);
  PetscViewerDestroy(viewer);
  PetscPrintf(PETSC_COMM_WORLD,"done.\n");

    }



  VecDestroy(gH);
  VecDestroy(gO);
  VecDestroy(gHO);
  VecDestroy(dgH);
  VecDestroy(dgO);
  VecDestroy(dgHO);
  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(f);

  VecDestroy(dgOH);
  VecDestroy(gOH);
  VecDestroy(tH);
  VecDestroy(tO);
  VecDestroy(tHO);
  VecDestroy(tOH);

  // ExtractAxis(BHD, g, 0);


  BGY3dH2OData_free(BHD);


  return PETSC_NULL;
}


/* solve with product ansatz g=g0*dg */
Vec BGY3d_solve_3site(PData PD, Vec g_ini, int vdim)
{
  BGY3dH2OData BHD;
  Vec g0H, g0O, g0HO, dgH, dgO, dgHO, dg_new, dg_new2, f, gH, gO, gHO;
  Vec dgOH, gOH;
  Vec tH, tO, tHO, tOH;
  real a=0.9, damp, damp_LJ, damp_start=0.0;
  real aH, aHO, aO, a0=0.9, a1=0.5, count=0.0, zpad; //, norm;
  int max_iter=25, iter, mycount=0, upwards=0, namecount=0; // in_iter
  char nameO[20], nameH[20], nameHO[20];
  PetscScalar dgH_norm, dgO_norm, dgHO_norm, norm_tol=1.0e-2;
  PetscScalar dgH_old, dgHO_old, dgO_old;

  // PetscScalar dgOH_norm;
  PetscTruth kflg, load_flag;
  PetscViewer viewer;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Fourier ansatz...\n");

  PetscOptionsHasName(PETSC_NULL,"-pair",&kflg);
  if(kflg)
    {

      BHD = BGY3dH2OData_Pair_malloc(PD);
      BHD->rho_H = 2.*BHD->rho_H;
      if( r_HH <0)
	{
	  PetscPrintf(PETSC_COMM_WORLD,"Solvent not a 3-Site model!\n");
	  exit(1);
	}
    }
  else
    {
      exit(1);
      //BHD = BGY3dFourierData_malloc(PD);
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }

  /* read BGY3d specific things from command line */
  /* Mixing parameter */
  PetscOptionsGetReal(PETSC_NULL,"-lambda",&a0, PETSC_NULL);
  if(a>1 || a<0)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lambda out of range: lambda=%f\n",a);
      exit(1);
    }
   /* Get damp_start from command line*/
  PetscOptionsGetReal(PETSC_NULL,"-damp_start",&damp_start, PETSC_NULL);
  /* Number of total iterations */
  PetscOptionsGetInt(PETSC_NULL,"-max_iter",&max_iter, PETSC_NULL);
  /* norm_tol for convergence test */
  PetscOptionsGetReal(PETSC_NULL,"-norm_tol",&norm_tol, PETSC_NULL);
  /* Zeropad */
  PetscOptionsGetReal(PETSC_NULL,"-zpad",&zpad, PETSC_NULL);
  /*********************************/

  PetscPrintf(PETSC_COMM_WORLD,"lambda = %f\n",a);
  PetscPrintf(PETSC_COMM_WORLD,"tolerance = %e\n",norm_tol);
  PetscPrintf(PETSC_COMM_WORLD,"zpad = %f\n",zpad);
  PetscPrintf(PETSC_COMM_WORLD,"max_iter = %d\n",max_iter);

  DACreateGlobalVector(BHD->da, &gH);
  DACreateGlobalVector(BHD->da, &gO);
  DACreateGlobalVector(BHD->da, &gHO);
  DACreateGlobalVector(BHD->da, &dgH);
  DACreateGlobalVector(BHD->da, &dgO);
  DACreateGlobalVector(BHD->da, &dgHO);
  DACreateGlobalVector(BHD->da, &dg_new);
  DACreateGlobalVector(BHD->da, &dg_new2);
  DACreateGlobalVector(BHD->da, &f);

  DACreateGlobalVector(BHD->da, &dgOH);
  DACreateGlobalVector(BHD->da, &gOH);

  DACreateGlobalVector(BHD->da, &tH);
  DACreateGlobalVector(BHD->da, &tO);
  DACreateGlobalVector(BHD->da, &tHO);
  DACreateGlobalVector(BHD->da, &tOH);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix */
  InitializeLaplaceMatrix(BHD, zpad);
  /* Create KSP environment */
  InitializeKSPSolver(BHD);
#endif

  g0H=BHD->gH_ini;
  g0O=BHD->gO_ini;
  g0HO=BHD->gHO_ini;

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dgHO,0);
  VecSet(dgOH,0.0);


/*   VecSetRandom_H2O(dgH,1.0); */
/*   VecSetRandom_H2O(dgO,1.0); */
/*   VecSetRandom_H2O(dgHO,1.0); */


  /* load initial configuration from file ??? */
  PetscOptionsHasName(PETSC_NULL,"-loadH2O",&load_flag);
  if(load_flag)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      /* dgH */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgH.bin",
			    FILE_MODE_READ,&viewer);
      VecLoad(viewer,VECMPI, &dgH);
      PetscViewerDestroy(viewer);
      /* dgO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgO.bin",
			    FILE_MODE_READ,&viewer);
      VecLoad(viewer,VECMPI, &dgO);
      PetscViewerDestroy(viewer);
      /* dgHO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgHO.bin",
			    FILE_MODE_READ,&viewer);
      VecLoad(viewer,VECMPI, &dgHO);
      PetscViewerDestroy(viewer);
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
    }

  VecSet(dg_new,0.0);






  //ComputeH2O_g( gOH, g0HO, dgOH);


  for( damp=damp_start; damp <=damp_start; damp+=0.1)
    {
      if(damp==-0.01)
	{
	  damp_LJ=0;
	  //a0=0.4;
	  RecomputeInitialData(BHD, 0, 1.0);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}
      else if(damp==0.0)
	{
	  damp_LJ=1.0;
	  //a0=0.5;
	  RecomputeInitialData(BHD, 0, 1.0);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}
      else
	{
	  damp_LJ=1.0;
	  //a0=0.0002/damp;
	  count+=1.0;
	  //a0 = 0.1/4./count;
	  //a0 = 0.1/count;
	  RecomputeInitialData(BHD, (damp), 1.0);
	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);
	}

/*       damp=0.05; */
/*       RecomputeInitialData(BHD, damp, 1.0); */

/*       Smooth_Function(BHD, g0HO, SL, SR, 0.0); */
/*       Smooth_Function(BHD, g0O, SL, SR, 0.0); */
/*       Smooth_Function(BHD, g0H, SL, SR, 0.0); */
/*       Zeropad_Function(BHD, g0HO, zpad, 0.0); */
/*       Zeropad_Function(BHD, g0O, zpad, 0.0); */
/*       Zeropad_Function(BHD, g0H, zpad, 0.0); */
      ImposeLaplaceBoundary(BHD, g0H, tH, BHD->xH, zpad, NULL);
      ImposeLaplaceBoundary(BHD, g0O, tH, BHD->xO, zpad, NULL);
      ImposeLaplaceBoundary(BHD, g0HO, tH, BHD->xHO, zpad, NULL);
      Zeropad_Function(BHD, g0HO, zpad, 0.0);
      Zeropad_Function(BHD, g0O, zpad, 0.0);
      Zeropad_Function(BHD, g0H, zpad, 0.0);
      /* g=g0*exp(-dg) */

      ComputeH2O_g( gHO, g0HO, dgHO);
      ComputeH2O_g( gH, g0H, dgH);
      ComputeH2O_g( gO, g0O, dgO);

      a1=a0;
      a=a0;
      aH=a;
      aO=a;
      aHO=a;

  for(iter=0; iter<max_iter; iter++)
    {

/*       Compute_dg_H2O_inter(BHD,  */
/* 			   BHD->fHO, BHD->fHO_l, gHO, gHO,  */
/* 			   BHD->ucHO_fft, BHD->rho_H,  */
/* 			   BHD->fO, BHD->fO_l, gO, gO,  */
/* 			   BHD->ucO_fft, BHD->rho_O,  */
/* 			   dg_new, f); */
/*       Compute_dg_H2O_inter(BHD, BHD->fH, BHD->fH_l, gH, gHO,  */
/* 			   BHD->ucH_fft, BHD->rho_H,  */
/* 			   BHD->fHO, BHD->fHO_l, gHO, gO,  */
/* 			   BHD->ucHO_fft, BHD->rho_O,  */
/* 			   dg_new, f); */

      if( !(iter%10) && iter>0 )
	{
	  a=a1;
	  //PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a);
	  aH=a;
	  aO=a;
	  aHO=a;
	}
      else if( iter==20)
	{
/*  	  a=0.1;  */
/*  	  PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a);  */
	  aH=a;
	  aO=a;
	  aHO=a;
	}
      else
	{
	  //a=0.01;
	  //PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a);
	  a=a0;
	  aH=a;
	  aO=a;
	  aHO=a;
	}

      //PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms: %e %e ", iter+1, NORM_REG, NORM_REG2);
      PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms:\t", iter+1);
      /* f=integral(g) */
      if(kflg)
	{
/* 	  Compute_dg_H2O_normalization_inter( BHD, gH, gH, tHO, f); */

/* 	  VecView(tHO,PETSC_VIEWER_STDERR_WORLD);          */
/*  	  exit(1);     */

	  //if(iter < 50) goto gH;
/* 	  if( !(iter%50) && iter>0 && NORM_REG>=2.0e-2) */
/* 	    NORM_REG/=2.; */

/* 	  ComputeH2O_g( gH, g0H, dgH); */
/* 	  ComputeH2O_g( gO, g0O, dgO); */
/* 	  ComputeH2O_g( gHO, g0HO, dgHO); */



	  goto gOH;
	  /* g_HO */
	  Compute_dg_H2O_inter(BHD,
			       BHD->fHO, BHD->fHO_l, gHO, gH,
			       BHD->ucHO_fft, BHD->rho_O, BHD->ucHO_0,
			       BHD->fH, BHD->fH_l, gH, gHO,
			       BHD->ucH_fft, BHD->rho_H, BHD->ucH_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);
	  //VecSet(dg_new,0.0);

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);          */
/*  	  exit(1);     */

	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gO, tO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gO, tO, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);
	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HH, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tHO, r_HH, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gHO, tHO, r_HH, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gHO, tHO, r_HH, dg_new2, f);
	  VecAXPY(dg_new, 1.0, dg_new2);



#ifdef INTRA1
	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fH, BHD->fH_l, tH, PETSC_NULL,
			       BHD->ucH_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gH, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 2.0, dg_new2);
#endif
#ifdef INTRA2
	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gH, r_HO, tHO, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fH, BHD->fH_l, tH, tHO,
				 BHD->ucH_fft, r_HO, dg_new2, f);
	  VecAXPY(dg_new, 2.0, dg_new2);
#endif

/* 	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);         */
/* 	  exit(1);    */

	  VecAXPY(dg_new, PD->beta, BHD->ucHO);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);
	  if(iter>=0)
	    {
	      ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xHO, zpad, NULL);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	    }

/* 	  VecNorm(dg_new, NORM_2, &norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"nrom=%e  ",norm); */

	  VecCopy(dgHO, f);
	  //VecAXPBY(dgHO, a, (1-a), dg_new);
	  VecAXPBY(dgHO, aHO, (1-aHO), dg_new);
	  VecAXPY(f, -1.0, dgHO);
	  VecNorm(f, NORM_INFINITY, &dgHO_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"HO= %e  (%f)  ",  dgHO_norm/aHO, aHO);
	  ComputeH2O_g( gHO, g0HO, dgHO);

/* 	  for(in_iter=0; in_iter<0; in_iter++) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD,"in_iter %d= ",in_iter); */

/* 	  Solve_NormalizationH2O(BHD, gH,  gO, gHO,  gOH, */
/* 				 tH, tO, tHO, tOH, dg_new2, f); */
	  goto gH;
	gOH:
	  //goto gH;
	  /* g_OH */
	  Compute_dg_H2O_inter(BHD,
			       BHD->fO, BHD->fO_l, gO, gHO,
			       BHD->ucO_fft, BHD->rho_O, BHD->ucO_0,
			       BHD->fHO, BHD->fHO_l, gHO, gH,
			       BHD->ucHO_fft, BHD->rho_H, BHD->ucHO_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);          */
/*  	  exit(1);     */


	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gH, tH, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gH, tH, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gH, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
 	  VecAXPY(dg_new, 2.0, dg_new2);



#ifdef INTRA1

	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HH, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fHO, BHD->fHO_l, tHO, PETSC_NULL,
			       BHD->ucHO_fft, r_HH, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gHO, r_HH, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);

	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fO, BHD->fO_l, tO, PETSC_NULL,
			       BHD->ucO_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);


#endif
#ifdef INTRA2
	  /* tO = gHO/int(gHO wHH) */
	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HH, gHO, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HH, tHO, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fHO, BHD->fHO_l, tO, tHO,
				 BHD->ucHO_fft, r_HH, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HH, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);

	  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gO, r_HO, tHO, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fO, BHD->fO_l, tO, tHO,
				 BHD->ucO_fft, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);
#endif

/* 	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);         */
/* 	  exit(1);   */

	  VecAXPY(dg_new, PD->beta, BHD->ucHO);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);

	  if(iter>=0)
	    {
	      ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xHO, zpad, NULL);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	    }
/* 	  VecNorm(dg_new, NORM_2, &norm); */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"nrom=%e  ",norm); */

	  VecCopy(dgHO, f);
	  //VecAXPBY(dgHO, a, (1-a), dg_new);
	  VecAXPBY(dgHO, aHO, (1-aHO), dg_new);
	  VecAXPY(f, -1.0, dgHO);
	  VecNorm(f, NORM_INFINITY, &dgHO_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"HO= %e  (%f)  ",  dgHO_norm/aHO, aHO);
	  //ComputeH2O_g( gHO, g0HO, dgHO);

	  //VecView(dgHO,PETSC_VIEWER_STDERR_WORLD);
 	  //exit(1);

	  //if(iter==max_iter-1) VecCopy(f, dgHO);
/* 	  } */
 	  /*****************************************/
/* 	  Solve_NormalizationH2O( BHD,  gH,  gO,  gHO, gOH, tH, tO, tHO, tOH, dg_new, f,  */
/* 			   norm_tol); */
	  /*****************************************/
/* 	  for(in_iter=0; in_iter<0; in_iter++) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD,"in_iter %d= ",in_iter); */
	gH:
	  //goto ende;
	  //if(damp==0.0 && iter<50) goto ende;
	  /* g_H */
  	  Compute_dg_H2O_inter(BHD,
			       BHD->fHO, BHD->fHO_l, gHO, gHO,
			       BHD->ucHO_fft, BHD->rho_O, BHD->ucHO_0,
			       BHD->fH, BHD->fH_l, gH, gH,
			       BHD->ucH_fft, BHD->rho_H, BHD->ucH_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);
	  //VecScale(dg_new, 0.5);

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);        */
/* 	  exit(1);   */

	  Solve_NormalizationH2O_smallII( BHD, gH, r_HH, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tH, r_HH, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gH, tH, r_HH, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gH, r_HH, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);

	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gHO, tHO, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);




#ifdef INTRA1

	  Solve_NormalizationH2O_smallII( BHD, gH, r_HH, gH, tH , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fH, BHD->fH_l, tH, PETSC_NULL,
			       BHD->ucH_fft, r_HH, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gH, r_HH, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);
 	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fHO, BHD->fHO_l, tHO, PETSC_NULL,
			       BHD->ucHO_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 1.0, dg_new2);

#endif
#ifdef INTRA2
	  /* tO = gH/int(gH wHH) */
	  Solve_NormalizationH2O_smallII( BHD, gH, r_HH, gH, tO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gH, r_HH, tH, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fH, BHD->fH_l, tO, tH,
				  BHD->ucH_fft, r_HH, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gH, r_HH, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);

	  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HO, tH, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fHO, BHD->fHO_l, tHO, tH,
				 BHD->ucHO_fft, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 1.0, dg_new2);
#endif

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);        */
/*   	  exit(1);   */


	  VecAXPY(dg_new, PD->beta, BHD->ucH);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);


	  if(iter>=0)
	    {
	      ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xH, zpad, NULL);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	    }
	  VecCopy(dgH, f);
	  //VecAXPBY(dgH, a, (1-a), dg_new);
	  VecAXPBY(dgH, aH, (1-aH), dg_new);
	  VecAXPY(f, -1.0, dgH);
	  VecNorm(f, NORM_INFINITY, &dgH_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"H= %e  (%f)  ", dgH_norm/aH, aH);
	  //ComputeH2O_g( gH, g0H, dgH);

	  //if(iter<50 )goto ende;
	  //if(iter==max_iter-1) VecCopy(f, dgH);
	  //ComputeH2O_Renormalization(BHD, gH);
/* 	  } */
/* 	  for(in_iter=0; in_iter<300; in_iter++) { */
/* 	    PetscPrintf(PETSC_COMM_WORLD,"in_iter %d= ",in_iter); */
	gO:
	  /* g_O */
	  //goto ende;
	  Compute_dg_H2O_inter(BHD,
			       BHD->fHO, BHD->fHO_l, gHO, gHO,
			       BHD->ucHO_fft, BHD->rho_H, BHD->ucHO_0,
			       BHD->fO, BHD->fO_l, gO, gO,
			       BHD->ucO_fft, BHD->rho_O, BHD->ucO_0,
			       dg_new, f);
	  VecScale(dg_new,damp_LJ);

/*  	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);         */
/*  	  exit(1);    */

	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnII(BHD, gHO, tHO, r_HO, dg_new2, f);
	  //Compute_dg_H2O_intra_lnIII(BHD, gHO, tHO, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  //Solve_NormalizationH2O_small( BHD, gO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 2.0, dg_new2);

/* 	  VecView(dg_new2,PETSC_VIEWER_STDERR_WORLD);        */
/* 	  exit(1);   */


#ifdef INTRA1

	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_intra(BHD, BHD->fHO, BHD->fHO_l, tHO, PETSC_NULL,
			       BHD->ucHO_fft, r_HO, dg_new2, f);
	  Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);
	  VecAXPY(dg_new, 2.0, dg_new2);

#endif
#ifdef INTRA2

	  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , dg_new2, f, zpad);
	  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HO, tO, f);
	  Compute_dg_H2O_intraIII(BHD, BHD->fHO, BHD->fHO_l, tHO, tO,
				 BHD->ucHO_fft, r_HO, dg_new2, f);
	  //Solve_NormalizationH2O_small( BHD, gHO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  //Solve_NormalizationH2O_smallII( BHD, gO, r_HO, dg_new2, dg_new2 , tOH, f, zpad);//
	  VecAXPY(dg_new, 2.0, dg_new2);
#endif

/* 	  VecView(dg_new,PETSC_VIEWER_STDERR_WORLD);       */
/* 	  exit(1);  */

	  VecAXPY(dg_new, PD->beta, BHD->ucO);
	  //Smooth_Function(BHD, dg_new, SL, SR, 0.0);
	  if(iter>=0)
	    {
	      ImposeLaplaceBoundary(BHD, dg_new, tH, BHD->xO, zpad, NULL);
	      Zeropad_Function(BHD, dg_new, zpad, 0.0);
	    }
	  VecCopy(dgO, f);
	  //VecAXPBY(dgO, a, (1-a), dg_new);
 	  VecAXPBY(dgO, aO, (1-aO), dg_new);
	  VecAXPY(f, -1.0,  dgO);
	  VecNorm(f, NORM_INFINITY, &dgO_norm);
	  PetscPrintf(PETSC_COMM_WORLD,"O= %e  (%f)  ", dgO_norm/aO, aO);
	  //ComputeH2O_g( gO, g0O, dgO);

	ende:
	  ;
	  ComputeH2O_g( gHO, g0HO, dgHO);
	  ComputeH2O_g( gH, g0H, dgH);
	  ComputeH2O_g( gO, g0O, dgO);

/* 	  CheckMax(gHO, "HO", 4.0); */
/* 	  CheckMax(gH, "H", 4.0); */
/* 	  CheckMax(gO, "O", 4.0); */

/* 	  norm = ComputeCharge(BHD, gH, gHO); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " H-OH=%e ", norm); */
/* 	  norm = ComputeCharge(BHD, gO, gHO); */
/* 	  PetscPrintf(PETSC_COMM_WORLD, " O-OH=%e ", norm); */


	  //if(iter==max_iter-1) VecCopy(f, dgO);
	  //ComputeH2O_Renormalization(BHD, gO);
/* 	  } */

/* 	  real max; */
/* 	  PetscScalar *f_vec; */
/* 	  VecGetArray(f, &f_vec); */
/* 	  max =0; */
/* 	  for(in_iter=0; in_iter<PD->N[0]*PD->N[1]*PD->N[2]; in_iter++) */
/* 	    { */
/* 	      if( fabs(f_vec[in_iter]) > max) */
/* 		{ */
/* 		  max = fabs(f_vec[in_iter]); */
/* 		  max_iter=in_iter; */
/* 		} */
/* 	    } */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"max at %d :%e %e ||",max_iter, max, f_vec[max_iter]); */
/* 	  VecRestoreArray(f, &f_vec); */



	}
      else {
        // nothing
      }

      //PetscPrintf(PETSC_COMM_WORLD,"\n");

      /* (fancy) step size control */
      mycount++;
      if( ((iter-1)%10) &&
	  (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
	{
	  upwards=1;
	}
      else if(iter>20 && !((iter-1)%10) && upwards==0 &&
	      (dgH_old<dgH_norm/aH || dgO_old<dgO_norm/aO || dgHO_old<dgHO_norm/aHO) )
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
      dgH_old = dgH_norm/aH;
      dgO_old = dgO_norm/aO;
      dgHO_old = dgHO_norm/aHO;
      /*********************************/

       //a = a0 * ( 1.0 + 0.5 * (0.5-(rand()/(real)RAND_MAX)));

      /* Set a */
/*       dgHO_norm /= a; */
/*       dgH_norm /= a; */
/*       dgO_norm /= a; */
/*       if(iter>0) */
/* 	{  */
/* 	  if( dgH_norm/dgH_old > ITERTOL_DOWN && dgO_norm/dgO_old > ITERTOL_DOWN && */
/* 	      dgHO_norm/dgHO_old > ITERTOL_DOWN && a<0.5) */
/* 	    a *= 1.1; */
/* 	  if ( (dgH_norm > dgH_old || dgO_norm > dgO_old || dgHO_norm > dgHO_old ) && a>0.01) */
/* 	    a /=2; */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"\tNew lambda= %f", a);  */

/* 	} */
/*       dgHO_old = dgHO_norm; */
/*       dgH_old = dgH_norm; */
/*       dgO_old = dgO_norm; */

      /* Set as */
/*       dgHO_norm /= aHO; */
/*       dgH_norm /= aH; */
/*       dgO_norm /= aO; */
/*       if(iter>0) */
/* 	{  */
/* 	  if( dgH_norm/dgH_old > ITERTOL_DOWN && dgH_norm/dgH_old <1 && aH<0.5 ) */
/* 	    aH *= 1.1; */
/* 	  if (dgH_norm > dgH_old&& aH>a ) */
/* 	    aH /= 2; */
/* 	  if(  dgO_norm/dgO_old > ITERTOL_DOWN && dgO_norm/dgO_old <1 && aO<0.5) */
/* 	    aO *= 1.1; */
/* 	  if (dgO_norm > dgO_old && aO>a) */
/* 	    aO /= 2; */
/* 	  if(  dgHO_norm/dgHO_old > ITERTOL_DOWN && dgHO_norm/dgHO_old <1 && aHO<0.5) */
/* 	    aHO *= 1.1; */
/* 	  if (dgHO_norm > dgHO_old && aHO>a) */
/* 	    aHO /= 2;      */
/* 	  PetscPrintf(PETSC_COMM_WORLD,"aHO= %f \taH= %f \taO= %f\n", aHO, aH, aO);  */
/* 	} */
/*       dgHO_old = dgHO_norm; */
/*       dgH_old = dgH_norm; */
/*       dgO_old = dgO_norm; */

/*       if(dgH_norm/a<=norm_tol &&  dgO_norm/a<=norm_tol && dgHO_norm/a<=norm_tol &&  */
/* 	 dgOH_norm/a<=norm_tol) */
/* 	break; */
      if(dgH_norm/aH<=norm_tol &&  dgO_norm/aO<=norm_tol && dgHO_norm/aHO<=norm_tol )//&&NORM_REG<=1.0e-2)
	break;

/*       VecSum(dgHO, &dgHO_norm); */
/*       VecSum(dgO, &dgO_norm); */
/*       VecSum(dgH, &dgH_norm); */
/*       PetscPrintf(PETSC_COMM_WORLD,"\t norms: dgHO= %e, dgH= %e, dgO= %e", dgHO_norm, dgH_norm, dgO_norm); */
      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }


/*   VecView(dg,PETSC_VIEWER_STDERR_WORLD); */

  /*************************************/
  /* output */
  namecount++;
  sprintf(nameH, "vecH-%d.m", namecount-1);
  sprintf(nameO, "vecO-%d.m", namecount-1);
  sprintf(nameHO, "vecHO-%d.m", namecount-1);
  PetscPrintf(PETSC_COMM_WORLD,"Writing files...");
  /* g_H */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameH,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecH.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gH,viewer);
  PetscViewerDestroy(viewer);
  /* g_b */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameO,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecO.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gO,viewer);
  PetscViewerDestroy(viewer);
  /* g_HO */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameHO,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gHO,viewer);
  PetscViewerDestroy(viewer);
  /* g_ba */
  //PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecba.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecab.m",FILE_MODE_WRITE,&viewer);
  //PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  //VecView(gba,viewer);
  //PetscViewerDestroy(viewer);
  PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  /************************************/
  /* save g to binary file */
  PetscOptionsHasName(PETSC_NULL,"-saveH2O",&load_flag);
  if(load_flag)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
      /* dgH */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgH.bin",
			    FILE_MODE_WRITE,&viewer);
      VecView(dgH,viewer);
      PetscViewerDestroy(viewer);
      /* dgO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgO.bin",
			    FILE_MODE_WRITE,&viewer);
      VecView(dgO,viewer);
      PetscViewerDestroy(viewer);
      /* dgHO */
      PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dgHO.bin",
			    FILE_MODE_WRITE,&viewer);
      VecView(dgHO,viewer);
      PetscViewerDestroy(viewer);
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
    }

  /************************************/
  /* save g2 to binary file */
  PetscPrintf(PETSC_COMM_WORLD,"Writing g2 files...");
  /* g2H */
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2H.bin",
			FILE_MODE_WRITE,&viewer);
  VecView(gH,viewer);
  PetscViewerDestroy(viewer);
  /* g2O */
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2O.bin",
			FILE_MODE_WRITE,&viewer);
  VecView(gO,viewer);
  PetscViewerDestroy(viewer);
  /* g2HO */
  PetscViewerBinaryOpen(PETSC_COMM_WORLD,"g2HO.bin",
			FILE_MODE_WRITE,&viewer);
  VecView(gHO,viewer);
  PetscViewerDestroy(viewer);
  PetscPrintf(PETSC_COMM_WORLD,"done.\n");

    }



  VecDestroy(gH);
  VecDestroy(gO);
  VecDestroy(gHO);
  VecDestroy(dgH);
  VecDestroy(dgO);
  VecDestroy(dgHO);
  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(f);

  VecDestroy(dgOH);
  VecDestroy(gOH);
  VecDestroy(tH);
  VecDestroy(tO);
  VecDestroy(tHO);
  VecDestroy(tOH);

  // ExtractAxis(BHD, g, 0);


  BGY3dH2OData_free(BHD);


  return PETSC_NULL;
}

