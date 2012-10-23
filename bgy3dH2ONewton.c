/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2ONewton.c,v 1.14 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"
#include "bgy3dH2OS.h"          /* for this */
#include "bgy3dH2O.h"
#include "bgy3dH2ONewton.h"

typedef struct H2Odg
{
  PetscScalar dgH, dgO, dgHO;
} H2Odg;

/*===================================================*/
/* Parametrisierung TIP3P :  */
/* sigma_H = 0.400   epsilon_H = 0.046  q_H = 0.417*/
/* sigma_O = 3.1506  epsilon_O = 0.1521 q_O =-0.834 */
/* r_OH = 0.9572 */
/* r_HH = 1.5139 */
/* theta_HOH = 104.52 */
/* mass H2O = 18.0154 u */
/* density = 1 kg/l = 0.6022142 u/A^3 => 0.033427745 / A^3 */
/* temperature : T= 298,15 K (25 C) => 0.5921 , => beta =1.6889 */

/* EPSILON0INV */
/* You have: e^2/4/pi/epsilon0/angstrom */
/* You want: kcal/avogadro/mol */
/* => 331.84164 */

#define sH 0.4 //0.400
#define eH 0.15 //0.046 //0.046
#define qH 0.417
#define sO 3.1506 //2.7165 //3.1506
#define eO 0.15 //0.1521 //0.1521
#define qO -0.834

#define r_HH  1.5139
#define r_HO  0.9572

#define EPSILON0INV 331.84164 //331.84164


static State *BGY3dH2OData_Pair_Newton_malloc(const ProblemData *PD)
{
  State *BHD;
  real interval[2], h[3], r[3], r_s, beta;
  int i[3], x[3], n[3];
  PetscScalar ***(fH_vec[3]),***(fO_vec[3]),***(fHO_vec[3]);
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]),***(fHOl_vec[3]);
  PetscScalar ***gHini_vec, ***gOini_vec, ***gHOini_vec;
  H2Odg ***pre_vec;
  real epsilonH, epsilonO, epsilonHO;
  real sigmaH, sigmaO, sigmaHO;
  real q2H, q2O, q2HO;

  BHD = (State*) malloc(sizeof(*BHD));

  /****************************************************/
  /* set Lennard-Jones and Coulomb parameters */
  /****************************************************/

  /* water hydrogen */
  BHD->LJ_paramsH[0] = eH;  /* epsilon  */
  BHD->LJ_paramsH[1] = sH;  /* sigma    */
  BHD->LJ_paramsH[2] = SQR(qH); /* charge product */
  epsilonH = BHD->LJ_paramsH[0];
  sigmaH = BHD->LJ_paramsH[1];
  q2H = BHD->LJ_paramsH[2];

  /* water oxygen */
  BHD->LJ_paramsO[0] = eO;  /* epsilon  */
  BHD->LJ_paramsO[1] = sO;  /* sigma    */
  BHD->LJ_paramsO[2] = SQR(qO); /* charge product */
  epsilonO = BHD->LJ_paramsO[0];
  sigmaO = BHD->LJ_paramsO[1];
  q2O = BHD->LJ_paramsO[2];

  /* water O-H mixed parameters */
  BHD->LJ_paramsHO[0] = sqrt(eH*eO);  /* epsilon  */
  BHD->LJ_paramsHO[1] = 0.5*(sH+sO);  /* sigma    */
  BHD->LJ_paramsHO[2] = qH*qO; /* charge product */
  epsilonHO = BHD->LJ_paramsHO[0];
  sigmaHO = BHD->LJ_paramsHO[1];
  q2HO = BHD->LJ_paramsHO[2];

  /****************************************************/

  BHD->PD = PD;

  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", PD->beta);
  /******************************/
  BHD->rhos[0] = 2.*PD->rho;
  BHD->rhos[1] = PD->rho;
  beta = PD->beta;

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  FOR_DIM
    h[dim]=PD->h[dim];

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &BHD->fft_mat, &BHD->da, &BHD->dc);
  const DA da = BHD->da;

#ifdef L_BOUNDARY
  /* Create Matrix with appropriate non-zero structure */
  DAGetMatrix (da, MATMPIAIJ, &BHD->M);
#endif

  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* Create  Petsc  Distributed  Array   with  3  degrees  of  freedom
     according to fftw data distribution: */
  {
    /* Get number of processes */
    int np;
    MPI_Comm_size (PETSC_COMM_WORLD, &np);

    int lx[1] = {PD->N[2]};
    int ly[1] = {PD->N[1]};
    int lz[np];

    assert (0);                 /* FIXME: untested, n[0]? */
    MPI_Allgather (&n[0], 1, MPI_INT, lz, 1, MPI_INT, PETSC_COMM_WORLD);

    {
      PetscInt dim, M, N, P, m, n, p, dof, s;
      DAPeriodicType wrap;
      DAStencilType st;
      DAGetInfo (da, &dim, &M, &N, &P, &m, &n, &p, &dof, &s, &wrap, &st);
      DACreate3d (PETSC_COMM_WORLD,
                  DA_NONPERIODIC,
                  DA_STENCIL_STAR ,
                  M, N, P,
                  m, n, p,
                  3, 0,
                  lx, ly, lz,
                  &BHD->da_newton);
    }
  }

  /* Create global vectors */
  DACreateGlobalVector(da, &(BHD->g_ini[0]));
  DACreateGlobalVector(da, &(BHD->g_ini[1]));
  DACreateGlobalVector(da, &(BHD->gHO_ini));
  DACreateGlobalVector(da, &(BHD->gH));
  DACreateGlobalVector(da, &(BHD->gO));
  DACreateGlobalVector(da, &(BHD->gHO));
  DACreateGlobalVector(da, &(BHD->dgH));
  DACreateGlobalVector(da, &(BHD->dgO));
  DACreateGlobalVector(da, &(BHD->dgHO));
  DACreateGlobalVector(da, &(BHD->u2[0][0]));
  DACreateGlobalVector(da, &(BHD->u2[1][1]));
  DACreateGlobalVector(da, &(BHD->u2[0][1]));
  DACreateGlobalVector(da, &(BHD->f));
  DACreateGlobalVector(da, &(BHD->f2));
  DACreateGlobalVector(da, &(BHD->f3));
  DACreateGlobalVector(da, &(BHD->f4));

  /* Preconditioner */
  DACreateGlobalVector(BHD->da_newton, &(BHD->pre));
  FOR_DIM
    {
      DACreateGlobalVector(da, &(BHD->F[0][0][dim]));
      DACreateGlobalVector(da, &(BHD->F[1][1][dim]));
      DACreateGlobalVector(da, &(BHD->F[0][1][dim]));
      DACreateGlobalVector(da, &(BHD->F_l[0][0][dim]));
      DACreateGlobalVector(da, &(BHD->F_l[1][1][dim]));
      DACreateGlobalVector(da, &(BHD->F_l[0][1][dim]));
      DACreateGlobalVector(da, &(BHD->v[dim]));
    }

  FOR_DIM
    {
      VecSet(BHD->F[0][0][dim],0.0);
      VecSet(BHD->F[1][1][dim],0.0);
      VecSet(BHD->F[0][1][dim],0.0);
      VecSet(BHD->F_l[0][0][dim],0.0);
      VecSet(BHD->F_l[1][1][dim],0.0);
      VecSet(BHD->F_l[0][1][dim],0.0);
    }
  VecSet(BHD->g_ini[0], 0.0);
  VecSet(BHD->g_ini[1], 0.0);
  VecSet(BHD->gHO_ini, 0.0);

  DAVecGetArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecGetArray(da, BHD->g_ini[1], &gOini_vec);
  DAVecGetArray(da, BHD->gHO_ini, &gHOini_vec);
  FOR_DIM
    {
      DAVecGetArray(da, BHD->F[0][0][dim], &(fH_vec[dim]));
      DAVecGetArray(da, BHD->F[1][1][dim], &(fO_vec[dim]));
      DAVecGetArray(da, BHD->F[0][1][dim], &(fHO_vec[dim]));
      DAVecGetArray(da, BHD->F_l[0][0][dim], &(fHl_vec[dim]));
      DAVecGetArray(da, BHD->F_l[1][1][dim], &(fOl_vec[dim]));
      DAVecGetArray(da, BHD->F_l[0][1][dim], &(fHOl_vec[dim]));
    }
  DAVecGetArray(BHD->da_newton, BHD->pre, (void*) &pre_vec);

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
	  gHini_vec[i[2]][i[1]][i[0]] +=
	    beta* Lennard_Jones( r_s, epsilonH, sigmaH);
	  gOini_vec[i[2]][i[1]][i[0]] +=
	    beta* Lennard_Jones( r_s, epsilonO, sigmaO);
	  gHOini_vec[i[2]][i[1]][i[0]] +=
	    beta* Lennard_Jones( r_s, epsilonHO, sigmaHO);

	  /* Coulomb short */
	  gHini_vec[i[2]][i[1]][i[0]] +=
	    beta* Coulomb_short( r_s, q2H);
	  gOini_vec[i[2]][i[1]][i[0]] +=
	    beta* Coulomb_short( r_s, q2O);
	  gHOini_vec[i[2]][i[1]][i[0]] +=
	    beta* Coulomb_short( r_s, q2HO);

	  /* Coulomb long */
/* 	  gHini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    beta* Coulomb_long( r_s, BHD->LJ_paramsH); */
/* 	  gOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    beta* Coulomb_long( r_s, BHD->LJ_paramsO); */
/* 	  gHOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */


	   /* Coulomb */
/* 	  gHini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    beta* Coulomb( r_s, BHD->LJ_paramsH); */
/* 	  gOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    beta* Coulomb( r_s, BHD->LJ_paramsO); */
/* 	  gHOini_vec[i[2]][i[1]][i[0]] +=  */
/* 	    beta* Coulomb( r_s, BHD->LJ_paramsHO); */

	   FOR_DIM
	    {
	      /* Lennard-Jones */
	      fH_vec[dim][i[2]][i[1]][i[0]] +=
		Lennard_Jones_grad( r_s, r[dim], epsilonH, sigmaH);
	      fO_vec[dim][i[2]][i[1]][i[0]] +=
		Lennard_Jones_grad( r_s, r[dim], epsilonO, sigmaO);
	      fHO_vec[dim][i[2]][i[1]][i[0]] +=
		Lennard_Jones_grad( r_s, r[dim], epsilonHO, sigmaHO);

 	      /* Coulomb short */
	      fH_vec[dim][i[2]][i[1]][i[0]] +=
		Coulomb_short_grad( r_s, r[dim], q2H);
	      fO_vec[dim][i[2]][i[1]][i[0]] +=
		Coulomb_short_grad( r_s, r[dim], q2O);
	      fHO_vec[dim][i[2]][i[1]][i[0]] +=
		Coulomb_short_grad( r_s, r[dim], q2HO);

	      /* Coulomb long */
/*  	      fHl_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fOl_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHOl_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Coulomb_long_grad( r_s, r[dim], BHD->LJ_paramsHO); */

	      /* Coulomb */
/* 	      fH_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Coulomb_grad( r_s, r[dim], BHD->LJ_paramsH); */
/* 	      fO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Coulomb_grad( r_s, r[dim], BHD->LJ_paramsO); */
/* 	      fHO_vec[dim][i[2]][i[1]][i[0]] +=  */
/* 		Coulomb_grad( r_s, r[dim], BHD->LJ_paramsHO); */


	    }

	   /* set Preconditioner */
	   pre_vec[i[2]][i[1]][i[0]].dgHO = exp(-gHOini_vec[i[2]][i[1]][i[0]]);
	   pre_vec[i[2]][i[1]][i[0]].dgH = exp(-gHini_vec[i[2]][i[1]][i[0]]);
	   pre_vec[i[2]][i[1]][i[0]].dgO = exp(-gOini_vec[i[2]][i[1]][i[0]]);

	}




  DAVecRestoreArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecRestoreArray(da, BHD->g_ini[1], &gOini_vec);
  DAVecRestoreArray(da, BHD->gHO_ini, &gHOini_vec);
  /* Preconditioner */
  DAVecRestoreArray(BHD->da_newton, BHD->pre, (void*) &pre_vec);
  FOR_DIM
    {
      DAVecRestoreArray(da, BHD->F[0][0][dim], &(fH_vec[dim]));
      DAVecRestoreArray(da, BHD->F[1][1][dim], &(fO_vec[dim]));
      DAVecRestoreArray(da, BHD->F[0][1][dim], &(fHO_vec[dim]));
      DAVecRestoreArray(da, BHD->F_l[0][0][dim], &(fHl_vec[dim]));
      DAVecRestoreArray(da, BHD->F_l[1][1][dim], &(fOl_vec[dim]));
      DAVecRestoreArray(da, BHD->F_l[0][1][dim], &(fHOl_vec[dim]));
    }

  /* Allocate memory for fft */
  DACreateGlobalVector (BHD->dc, &BHD->fft_scratch);
  /* FIXME: where is the cleanup code? */

  FOR_DIM
    DACreateGlobalVector (BHD->dc, &BHD->fg2_fft[dim]);

  DACreateGlobalVector (BHD->dc, &BHD->gfg2_fft);

  /* FIXME: these probably differ only by factors: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD->dc, &BHD->u2_fft[i][j]);
        BHD->u2_fft[j][i] = BHD->u2_fft[i][j];
      }

  /* Compute fft from Coulomb potential (long) */
  ComputeFFTfromCoulomb(BHD, BHD->u2[0][1], BHD->F_l[0][1], BHD->u2_fft[0][1],
			q2HO);
  ComputeFFTfromCoulomb(BHD, BHD->u2[0][0], BHD->F_l[0][0], BHD->u2_fft[0][0],
			q2H);
  ComputeFFTfromCoulomb(BHD, BHD->u2[1][1], BHD->F_l[1][1], BHD->u2_fft[1][1],
			q2O);

  FOR_DIM
    {
      VecAXPY(BHD->F[0][1][dim], 1.0, BHD->F_l[0][1][dim]);
      VecAXPY(BHD->F[0][0][dim], 1.0, BHD->F_l[0][0][dim]);
      VecAXPY(BHD->F[1][1][dim], 1.0, BHD->F_l[1][1][dim]);
	}

  return BHD;
}




// FIXME: warning of snes cannot be eliminated becaused this would be passed
//      to SNESSetFunction()
static PetscErrorCode ComputeH2OFunction(SNES snes, Vec u, Vec f, void *data)
{
  (void) snes;                  /* FIXME: interface obligation? */

  State *BHD;
  H2Odg ***dg_struct;
  PetscScalar ***dgH_vec, ***dgHO_vec, ***dgO_vec;
  int i[3], x[3], n[3];
  Vec dgH, dgHO, dgO, gH, gHO, gO, help;
  Vec tHO, tO, tH, help2;

  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "--- Function evaluation starts...\n");

  BHD = (State*) data;
  gH= BHD->gH;
  gO= BHD->gO;
  gHO= BHD->gHO;
  dgH= BHD->dgH;
  dgO= BHD->dgO;
  dgHO= BHD->dgHO;
  help = BHD->f2;
  help2= BHD->f3;
  tHO =  BHD->f4;
  tO  =  BHD->f4;
  tH  =  BHD->f4;
  const real zpad = BHD->PD->zpad;

  /* Get arrays from PETSC Vectors */
  DAVecGetArray(BHD->da, dgH, &dgH_vec);
  DAVecGetArray(BHD->da, dgHO, &dgHO_vec);
  DAVecGetArray(BHD->da, dgO, &dgO_vec);
  DAVecGetArray(BHD->da_newton, u, (void*) &dg_struct);

  DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from u to single Vectors */
	  dgHO_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgHO;
	  dgH_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgH;
	  dgO_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgO;
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray(BHD->da, dgH, &dgH_vec);
  DAVecRestoreArray(BHD->da, dgHO, &dgHO_vec);
  DAVecRestoreArray(BHD->da, dgO, &dgO_vec);
  DAVecRestoreArray(BHD->da_newton, u, (void*) &dg_struct);

  /* Compute g's from dg's */
  Zeropad_Function(BHD, dgO, zpad, 0.0);
  Zeropad_Function(BHD, dgH, zpad, 0.0);
  Zeropad_Function(BHD, dgHO, zpad, 0.0);
  ComputeH2O_g( gHO, BHD->gHO_ini, dgHO);
  ComputeH2O_g( gH,  BHD->g_ini[0] , dgH);
  ComputeH2O_g( gO,  BHD->g_ini[1] , dgO);


  /* Compute right hand side */
  /***********************************************************/
  /* gOH */
  /***********************************************************/
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "Computing dgHO, ");
  VecCopy(dgHO, help);
  Compute_dg_H2O_inter(BHD,
		       BHD->F[1][1], BHD->F_l[1][1], gO, gHO,
		       BHD->u2_fft[1][1], BHD->rhos[1],
		       BHD->F[0][1], BHD->F_l[0][1], gHO, gH,
		       BHD->u2_fft[0][1], BHD->rhos[0],
		       dgHO, BHD->f);
  VecAXPY(dgHO, BHD->PD->beta, BHD->u2[0][1]);
  /************************************************************/
  /* intra molecular part */
  /************************************************************/
  //goto gHO_end;
  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gH, tH , help2, BHD->f, zpad);
  Compute_dg_H2O_intra_ln(BHD, tH, r_HO, help2);
  VecCopy (help2, BHD->f); /* FIXME: need that? */
  VecAXPY(dgHO, 2.0, help2);
#ifdef INTRA1
  Solve_NormalizationH2O_smallII( BHD, gHO, r_HH, gHO, tHO , help2, BHD->f, zpad);
  Compute_dg_H2O_intra(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, PETSC_NULL,
		       BHD->u2_fft[0][1], r_HH, help2, BHD->f);
  Solve_NormalizationH2O_small( BHD, gHO, r_HH, help2, help2 , tHO, BHD->f, zpad);
  VecAXPY(dgHO, 1.0, help2);
  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , help2, BHD->f, zpad);
  Compute_dg_H2O_intra(BHD, BHD->F[1][1], BHD->F_l[1][1], tO, PETSC_NULL,
		       BHD->u2_fft[1][1], r_HO, help2, BHD->f);
  Solve_NormalizationH2O_small( BHD, gO, r_HO, help2, help2 , tHO, BHD->f, zpad);
  VecAXPY(dgHO, 1.0, help2);
#endif
#ifdef INTRA2
  /* tO = gHO/int(gHO wHH) */
  Solve_NormalizationH2O_smallII( BHD, gHO, r_HH, gHO, tO , help2, BHD->f, zpad);
  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HH, tHO, BHD->f);
  Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tO, tHO,
			 BHD->u2_fft[0][1], r_HH, help2, BHD->f);
  VecAXPY(dgHO, 1.0, help2);
  Compute_dg_H2O_normalization_intra( BHD, gO, r_HO, tHO, BHD->f);
  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gO, tO , help2, BHD->f, zpad);
  Compute_dg_H2O_intraIII(BHD, BHD->F[1][1], BHD->F_l[1][1], tO, tHO,
			 BHD->u2_fft[1][1], r_HO, help2, BHD->f);
  VecAXPY(dgHO, 1.0, help2);
#endif
  /***********************************************************/
  //gHO_end:
  ImposeLaplaceBoundary(BHD, dgHO,BHD->v[0], BHD->v[1], zpad, NULL);
  Zeropad_Function(BHD, dgHO, zpad, 0.0);
  VecAYPX(dgHO, -1.0, help);
  /***********************************************************/
  /* gH */
  /***********************************************************/
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "dgH, ");
  VecCopy(dgH, help);
  Compute_dg_H2O_inter(BHD,
		       BHD->F[0][1], BHD->F_l[0][1], gHO, gHO,
		       BHD->u2_fft[0][1], BHD->rhos[1],
		       BHD->F[0][0], BHD->F_l[0][0], gH, gH,
		       BHD->u2_fft[0][0], BHD->rhos[0],
		       dgH, BHD->f);
  VecAXPY(dgH, BHD->PD->beta, BHD->u2[0][0]);
  /************************************************************/
  /* intra molecular part */
  /************************************************************/
  //goto gH_end;
  Solve_NormalizationH2O_smallII( BHD, gH, r_HH, gH, tH , help2, BHD->f, zpad);
  Compute_dg_H2O_intra_ln(BHD, tH, r_HH, help2);
  VecCopy (help2, BHD->f); /* FIXME: need that? */
  VecAXPY(dgH, 1.0, help2);
  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , help2, BHD->f, zpad);
  Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, help2);
  VecCopy (help2, BHD->f); /* FIXME: need that? */
  VecAXPY(dgH, 1.0, help2);

#ifdef INTRA1
  Solve_NormalizationH2O_smallII( BHD, gH, r_HH, gH, tH , help2, BHD->f, zpad);
  Compute_dg_H2O_intra(BHD, BHD->F[0][0], BHD->F_l[0][0], tH, PETSC_NULL,
		       BHD->u2_fft[0][0], r_HH, help2, BHD->f);
  Solve_NormalizationH2O_small( BHD, gH, r_HH, help2, help2 , tH, BHD->f, zpad);
  VecAXPY(dgH, 1.0, help2);
  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , help2, BHD->f, zpad);
  Compute_dg_H2O_intra(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, PETSC_NULL,
		       BHD->u2_fft[0][1], r_HO, help2, BHD->f);
  Solve_NormalizationH2O_small( BHD, gHO, r_HO, help2, help2 , tH, BHD->f, zpad);
  VecAXPY(dgH, 1.0, help2);
#endif
#ifdef INTRA2
  /* tO = gH/int(gH wHH) */
  Solve_NormalizationH2O_smallII( BHD, gH, r_HH, gH, tO , help2, BHD->f, zpad);
  Compute_dg_H2O_normalization_intra( BHD, gH, r_HH, tH, BHD->f);
  Compute_dg_H2O_intraIII(BHD, BHD->F[0][0], BHD->F_l[0][0], tO, tH,
			 BHD->u2_fft[0][0], r_HH, help2, BHD->f);
  VecAXPY(dgH, 1.0, help2);
  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HO, tH, BHD->f);
  Solve_NormalizationH2O_smallII( BHD, gH, r_HO, gHO, tHO , help2, BHD->f, zpad);
  Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, tH,
			 BHD->u2_fft[0][1], r_HO, help2, BHD->f);
  VecAXPY(dgH, 1.0, help2);
#endif
  /***********************************************************/
  //gH_end:
  ImposeLaplaceBoundary(BHD, dgH,BHD->v[0], BHD->v[1], zpad, NULL);
  Zeropad_Function(BHD, dgH, zpad, 0.0);
  VecAYPX(dgH, -1.0, help);
  /***********************************************************/
  /* gO */
  /***********************************************************/
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "dgO... ");
  VecCopy(dgO, help);
  Compute_dg_H2O_inter(BHD,
		       BHD->F[0][1], BHD->F_l[0][1], gHO, gHO,
		       BHD->u2_fft[0][1], BHD->rhos[0],
		       BHD->F[1][1], BHD->F_l[1][1], gO, gO,
		       BHD->u2_fft[1][1], BHD->rhos[1],
		       dgO, BHD->f);
  VecAXPY(dgO, BHD->PD->beta, BHD->u2[1][1]);
  /************************************************************/
  /* intra molecular part */
  /************************************************************/
  //goto gO_end;
  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , help2, BHD->f, zpad);
  Compute_dg_H2O_intra_ln(BHD, tHO, r_HO, help2);
  VecCopy (help2, BHD->f); /* FIXME: need that? */
  VecAXPY(dgO, 2.0, help2);
#ifdef INTRA1
  Solve_NormalizationH2O_smallII( BHD, gHO, r_HO, gHO, tHO , help2, BHD->f, zpad);
  Compute_dg_H2O_intra(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, PETSC_NULL,
		       BHD->u2_fft[0][1], r_HO, help2, BHD->f);
  Solve_NormalizationH2O_small( BHD, gHO, r_HO, help2, help2 , tO, BHD->f, zpad);
  VecAXPY(dgO, 2.0, help2);
#endif
#ifdef INTRA2
  Solve_NormalizationH2O_smallII( BHD, gO, r_HO, gHO, tHO , help2, BHD->f, zpad);
  Compute_dg_H2O_normalization_intra( BHD, gHO, r_HO, tO, BHD->f);
  Compute_dg_H2O_intraIII(BHD, BHD->F[0][1], BHD->F_l[0][1], tHO, tO,
			 BHD->u2_fft[0][1], r_HO, help2, BHD->f);
  VecAXPY(dgO, 2.0, help2);
#endif
  /***********************************************************/
  //gO_end:
  ImposeLaplaceBoundary(BHD, dgO,BHD->v[0], BHD->v[1], zpad, NULL);
  Zeropad_Function(BHD, dgO, zpad, 0.0);
  VecAYPX(dgO, -1.0, help);
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, " done.\n");


/*   Smooth_Function(BHD, dgHO, SL, SR, 0.0); */
/*   Smooth_Function(BHD, dgH, SL, SR, 0.0); */
/*   Smooth_Function(BHD, dgO, SL, SR, 0.0); */


/*   ImposeLaplaceBoundary(BHD, dgH, BHD->v[0], BHD->v[1], zpad); */
/*   ImposeLaplaceBoundary(BHD, dgO, BHD->v[0], BHD->v[1], zpad); */
/*   ImposeLaplaceBoundary(BHD, dgHO,BHD->v[0], BHD->v[1], zpad); */
/*   Zeropad_Function(BHD, dgHO, zpad, 0.0); */
/*   Zeropad_Function(BHD, dgO, zpad, 0.0); */
/*   Zeropad_Function(BHD, dgH, zpad, 0.0); */

  /* Get arrays from PETSC Vectors */
  DAVecGetArray(BHD->da, dgH , &dgH_vec);
  DAVecGetArray(BHD->da, dgHO, &dgHO_vec);
  DAVecGetArray(BHD->da, dgO , &dgO_vec);
  DAVecGetArray(BHD->da_newton, f, (void*) &dg_struct);


  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from single Vectors to f */
	  dg_struct[i[2]][i[1]][i[0]].dgHO= dgHO_vec[i[2]][i[1]][i[0]];
	  dg_struct[i[2]][i[1]][i[0]].dgH= dgH_vec[i[2]][i[1]][i[0]];
	  dg_struct[i[2]][i[1]][i[0]].dgO= dgO_vec[i[2]][i[1]][i[0]];
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray(BHD->da, dgH , &dgH_vec);
  DAVecRestoreArray(BHD->da, dgHO, &dgHO_vec);
  DAVecRestoreArray(BHD->da, dgO , &dgO_vec);
  DAVecRestoreArray(BHD->da_newton, f, (void*) &dg_struct);

  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "--- Function evaluation finished.\n");


/*   WriteH2ONewtonPlain(BHD, f); */
/*   exit(1); */

  return 0;

}


static void WriteH2ONewtonSolution(State *BHD, Vec u)
{
  H2Odg ***dg_struct;
  PetscScalar ***dgH_vec, ***dgHO_vec, ***dgO_vec;
  int i[3], x[3], n[3];
  Vec dgH, dgHO, dgO, gH, gHO, gO;

  PetscPrintf(PETSC_COMM_WORLD,"Writing files...");

  gH= BHD->gH;
  gO= BHD->gO;
  gHO= BHD->gHO;
  dgH= BHD->dgH;
  dgO= BHD->dgO;
  dgHO= BHD->dgHO;


  /* Get arrays from PETSC Vectors */
  DAVecGetArray(BHD->da, dgH, &dgH_vec);
  DAVecGetArray(BHD->da, dgHO, &dgHO_vec);
  DAVecGetArray(BHD->da, dgO, &dgO_vec);
  DAVecGetArray(BHD->da_newton, u, (void*) &dg_struct);

  DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from u to single Vectors */
	  dgHO_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgHO;
	  dgH_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgH;
	  dgO_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgO;
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray(BHD->da, dgH, &dgH_vec);
  DAVecRestoreArray(BHD->da, dgHO, &dgHO_vec);
  DAVecRestoreArray(BHD->da, dgO, &dgO_vec);
  DAVecRestoreArray(BHD->da_newton, u, (void*) &dg_struct);

  /* Copmute g's from dg's */
  ComputeH2O_g( gHO, BHD->gHO_ini, dgHO);
  ComputeH2O_g( gH,  BHD->g_ini[0] , dgH);
  ComputeH2O_g( gO,  BHD->g_ini[1] , dgO);

  /*************************************/
  /* output */
  bgy3d_save_vec_ascii ("vec00.m", gH);
  bgy3d_save_vec_ascii ("vec11.m", gO);
  bgy3d_save_vec_ascii ("vec01.m", gHO);

  /* save g2 to binary file */
  PetscPrintf(PETSC_COMM_WORLD,"Writing g2 files...");
  bgy3d_save_vec ("g00.bin", gH);
  bgy3d_save_vec ("g11.bin", gO);
  bgy3d_save_vec ("g01.bin", gHO);
  PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  /************************************/

}



/* apply preconditioner matrix: diagonal scaling */
static PetscErrorCode ComputePreconditioner_H2O(void *data, Vec x, Vec y)
{
  State *BHD;
  PetscErrorCode ierr;

  BHD = (State*) data;
  ierr=VecPointwiseMult(y,BHD->pre,x);


/*   VecView(x,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

  return ierr;
}

/*#include "petscksp.h"
#include "/opt/packages/petsc/src/snes/snesimpl.h"
#include "/opt/packages/petsc/src/snes/mf/snesmfj.h"
#include "/opt/packages/petsc/src/mat/matimpl.h"
PetscErrorCode MatMult_MFFD(Mat mat,Vec a,Vec y);
PetscErrorCode PETSCSNES_DLLEXPORT MatCreate_MFFD(Mat A);
*/

Vec BGY3d_SolveNewton_H2O(const ProblemData *PD, Vec g_ini)
{
  SNES snes;
  KSP ksp;
  PC  pc;
  State *BHD;
  PetscTruth flg;
  Vec u, f, b, v1, v2;
  real damp;
  // Mat M;
  // int local_size;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Newton ...\n");

  BHD = BGY3dH2OData_Pair_Newton_malloc(PD);

  /* Get damp_start from command line*/
  const real damp_start = PD->damp;

  /* Zeropad */
  const real zpad = PD->zpad;

  DACreateGlobalVector(BHD->da_newton, &f);
  DACreateGlobalVector(BHD->da_newton, &u);
  DACreateGlobalVector(BHD->da_newton, &b);
  DACreateGlobalVector(BHD->da_newton, &v1);
  DACreateGlobalVector(BHD->da_newton, &v2);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix */
  InitializeLaplaceMatrix(BHD, zpad);
  /* Create KSP environment */
  InitializeKSPSolver(BHD);
#endif

  /* Create SNES environment */
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp, &pc);
  /* set rtol, atol, dtol, maxits */
  KSPSetTolerances(ksp, 1.0e-3, 1.0e-10, 1.0e+5, 1000);
  /* set atol, rtol, stol , its, fct. eval. */
  SNESSetTolerances(snes, 5.0e-2, 1.0e-5, 1.0e-4 , 50, 10000);
  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType(snes, SNESLS);
  flg = bgy3d_getopt_test ("--user-precond");
  if (flg) { /* user-defined precond */
    /* Set user defined preconditioner */
    PCSetType(pc,PCSHELL);
    PCShellSetApply(pc,ComputePreconditioner_H2O);
    PCShellSetContext(pc,BHD);
  } else
    /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
    PCSetType( pc, PCJACOBI);
  /* set function */
  SNESSetFunction(snes, f, ComputeH2OFunction, (void*)BHD);

  /* runtime options will override default parameters */
  SNESSetFromOptions(snes);

  /* set initial guess */
  VecSet(u, 0.0);
  //VecSetRandom_H2O(u, 0.5);

  //MatCreateSNESMF(snes, u, &M);


/*   VecGetLocalSize(u, &local_size); */
/*   MatCreateShell(PETSC_COMM_WORLD, local_size, local_size, 3*PD->N3, 3*PD->N3, (void*)BHD, &M); */
/*   MatCreate_MFFD(M); */
/*   MatShellSetContext(M, (void*)BHD); */
/*   MatShellSetOperation(M,MATOP_MULT,(void*)MatMult_MFFD); */
/*   KSPSetOperators(ksp, M, M, SAME_NONZERO_PATTERN); */
/*   KSPGetPC(ksp, &pc); */
/*   PCSetType( pc, PCNONE); */

/*   VecSet(b, 1.0); */
/*   KSPSolve(ksp, b, u); */
  //KSPInitialResidual(ksp, u, v1, v2, f, b);

  for(damp=damp_start; damp<=1; damp+= 0.01)
    {

      RecomputeInitialData(BHD, (damp), 1.0);
/*       Smooth_Function(BHD, BHD->gHO_ini, SL, SR, 0.0); */
/*       Smooth_Function(BHD, BHD->g_ini[0], SL, SR, 0.0); */
/*       Smooth_Function(BHD, BHD->g_ini[1], SL, SR, 0.0); */
/*       Zeropad_Function(BHD, BHD->gHO_ini, zpad, 0.0); */
/*       Zeropad_Function(BHD, BHD->g_ini[0], zpad, 0.0); */
/*       Zeropad_Function(BHD, BHD->g_ini[1], zpad, 0.0); */
      ImposeLaplaceBoundary(BHD, BHD->g_ini[0], BHD->v[0], BHD->v[1], zpad, NULL);
      ImposeLaplaceBoundary(BHD, BHD->g_ini[1], BHD->v[0], BHD->v[1], zpad, NULL);
      ImposeLaplaceBoundary(BHD, BHD->gHO_ini,BHD->v[0], BHD->v[1], zpad, NULL);
      Zeropad_Function(BHD, BHD->gHO_ini, zpad, 0.0);
      Zeropad_Function(BHD, BHD->g_ini[0], zpad, 0.0);
      Zeropad_Function(BHD, BHD->g_ini[1], zpad, 0.0);
      /* solve problem */
      SNESSolve(snes, PETSC_NULL, u);

      /* Get Solution */
      //SNESGetSolution(snes, &u);

      /* Write out solution */
      WriteH2ONewtonSolution(BHD, u);
    }




/*   VecSet(u, 0.0); */
/*   ComputeH2OFunction(snes, u, f, BHD); */
/*   WriteH2ONewtonSolution(BHD, f); */



  VecDestroy(f);
  VecDestroy(u);
  VecDestroy(b);
  SNESDestroy(snes);

  assert (0);                   /* FIXME: clean up BHD? */
  return PETSC_NULL;

}




