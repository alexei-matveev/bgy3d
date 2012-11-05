/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2OSNewton.c,v 1.4 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-getopt.h"
#include "bgy3d-impure.h"       /* for this */
#include "bgy3d-pure.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-poisson.h"      /* laplace staff */
#include "bgy3d-newton-impure.h"

typedef struct H2OSdg
{
  PetscScalar dgH, dgO;
} H2OSdg;

typedef struct H2OSdgF
{
  PetscScalar dgHre, dgHim, dgOre, dgOim;

} H2OSdgF;

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

#define sH 1.8 //0.4 //0.400
#define eH 0.15 //0.046 //0.046
#define qH 0.417
#define sO 1.8 //2.8509 //3.1506 //2.7165 //3.1506
#define eO 0.15 //0.1521 //0.1521
#define qO -0.834

#define r_HH  1.5139
#define r_HO  0.9572

static void RecomputeInitialSoluteData(State *BHD, real damp, real damp_LJ);

static State *BGY3dH2OData_Newton_malloc(const ProblemData *PD)
{
  State *BHD;
  int x[3], n[3];
  PetscErrorCode ierr;

  BHD = (State*) malloc(sizeof(*BHD));

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

  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  //PetscPrintf(PETSC_COMM_WORLD, "Boundary smoothing parameters : SL= %f  SR= %f\n", SL, SR);
  //PetscPrintf(PETSC_COMM_WORLD, "ZEROPAD= %f\n", ZEROPAD);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", PD->beta);
  /******************************/
  BHD->rhos[0] = 2.*PD->rho;
  BHD->rhos[1] = PD->rho;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &BHD->fft_mat, &BHD->da, &BHD->dc);
  const DA da = BHD->da;

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
    int err = MPI_Allgather (&n[0], 1, MPI_INT, lz, 1, MPI_INT, PETSC_COMM_WORLD);
    assert (!err);

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
                  2, 0,
                  lx, ly, lz,
                  &BHD->da_newton);
      DACreate3d (PETSC_COMM_WORLD,
                  DA_NONPERIODIC,
                  DA_STENCIL_STAR ,
                  M, N, P,
                  m, n, p,
                  4, 0,
                  lx, ly, lz,
                  &BHD->da_newtonF);
    }
  }

  /*
   * CHKERRQ()  is  a  macro that  expands  to  an  if with  a  return
   * statement  in  the if-block.   It  is  not  suitable for  use  in
   * functions  returning  anything   that  but  PetscErrorCode.  This
   * function returns State*.
   */

  /* Create global vectors */
  ierr = DACreateGlobalVector(da, &(BHD->g_ini[0])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g_ini[1])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->gHO_ini)); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->u2[0][0])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->u2[1][1])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->u2[0][1])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g2[0][0])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g2[1][1])); // CHKERRQ(ierr);
  assert (!ierr);
  ierr = DACreateGlobalVector(da, &(BHD->g2[0][1])); // CHKERRQ(ierr);
  assert (!ierr);
  DACreateGlobalVector(da, &(BHD->gH));
  DACreateGlobalVector(da, &(BHD->gO));
  DACreateGlobalVector(da, &(BHD->dgH));
  DACreateGlobalVector(da, &(BHD->dgO));
  DACreateGlobalVector(da, &(BHD->f));
  DACreateGlobalVector(da, &(BHD->f2));
  DACreateGlobalVector(da, &(BHD->f3));
  DACreateGlobalVector(da, &(BHD->f4));
    /* Preconditioner */
  DACreateGlobalVector(BHD->da_newton, &(BHD->pre));
  FOR_DIM
    {

      ierr = DACreateGlobalVector(da, &(BHD->F[0][0][dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->F[1][1][dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->F[0][1][dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->F_l[0][0][dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->F_l[1][1][dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->F_l[0][1][dim])); // CHKERRQ(ierr);
      assert (!ierr);
      ierr = DACreateGlobalVector(da, &(BHD->v[dim])); // CHKERRQ(ierr);
      assert (!ierr);
    }

  /* Allocate memory for fft */
  DACreateGlobalVector (BHD->dc, &BHD->fft_scratch);
  /* FIXME: where is the cleanup code? */

  FOR_DIM
    {
      DACreateGlobalVector (BHD->dc, &BHD->fs_g2_fft[0][0][dim]);
      DACreateGlobalVector (BHD->dc, &BHD->fs_g2_fft[1][1][dim]);
      DACreateGlobalVector (BHD->dc, &BHD->fs_g2_fft[0][1][dim]);
      DACreateGlobalVector (BHD->dc, &BHD->fg2_fft[dim]);
    }

  DACreateGlobalVector (BHD->dc, &BHD->gfg2_fft);

  /* FIXME: these probably differ only by factors: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD->dc, &BHD->u2_fft[i][j]);
        BHD->u2_fft[j][i] = BHD->u2_fft[i][j];
      }

  DACreateGlobalVector (BHD->dc, &BHD->wHO_fft);
  DACreateGlobalVector (BHD->dc, &BHD->wHH_fft);



  /* Read g^2  from file */
  ReadPairDistribution(BHD, "g2_OO", BHD->g2[1][1]);
  ReadPairDistribution(BHD, "g2_HH", BHD->g2[0][0]);
  ReadPairDistribution(BHD, "g2_HO", BHD->g2[0][1]);

  /* Compute initial data */
  //RecomputeInitialFFTs(BHD, 1.0, 1.0);

  /* Compute Solute dependent initial data */
  //RecomputeInitialSoluteData(BHD, 1.0, 1.0);

  return BHD;
}




static void InitializePreconditioner(State *BHD)
{
  DA da;
  H2OSdg ***pre_vec;
  PetscScalar ***gHini_vec, ***gOini_vec;
  int i[3], n[3], x[3];

  da = BHD->da;

  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  DAVecGetArray(BHD->da_newton, BHD->pre, (void*) &pre_vec);
  DAVecGetArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecGetArray(da, BHD->g_ini[1], &gOini_vec);

   /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{

	  /* set Preconditioner */
	  if( exp(-gHini_vec[i[2]][i[1]][i[0]]) <1.0e-6)
	    pre_vec[i[2]][i[1]][i[0]].dgH = 0;
	  else
	    pre_vec[i[2]][i[1]][i[0]].dgH = 1;
	   if( exp(-gOini_vec[i[2]][i[1]][i[0]]) <1.0e-6)
	    pre_vec[i[2]][i[1]][i[0]].dgO = 0;
	  else
	    pre_vec[i[2]][i[1]][i[0]].dgO = 1;

	}
  DAVecRestoreArray(BHD->da_newton, BHD->pre, (void*) &pre_vec);
  DAVecRestoreArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecRestoreArray(da, BHD->g_ini[1], &gOini_vec);

}


// FIXME: warning of snes cannot be eliminated becaused this would be passed
//      to SNESSetFunction()
static PetscErrorCode ComputeH2OSFunction(SNES snes, Vec u, Vec f, void *data)
{
  (void) snes;                  /* FIXME: interface obligation? */

  State *BHD;
  H2OSdg ***dg_struct;
  PetscScalar ***dgH_vec, ***dgO_vec;
  int i[3], x[3], n[3];
  Vec dgH, dgO, gH, gO, help;
  Vec tO, tH, help2, ff;

  static int counter=0;
  counter++;
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "--- Function evaluation starts...\n");

  BHD = (State*) data;
  gH= BHD->gH;
  gO= BHD->gO;
  dgH= BHD->dgH;
  dgO= BHD->dgO;
  help = BHD->f2;
  help2= BHD->f3;
  tO  =  BHD->f4;
  tH  =  BHD->f4;
  ff = BHD->f;

  /* Get arrays from PETSC Vectors */
  DAVecGetArray(BHD->da, dgH, &dgH_vec);
  DAVecGetArray(BHD->da, dgO, &dgO_vec);
  DAVecGetArray(BHD->da_newton, u, (void*) &dg_struct);

  DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from u to single Vectors */
	  dgH_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgH;
	  dgO_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgO;
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray(BHD->da, dgH, &dgH_vec);
  DAVecRestoreArray(BHD->da, dgO, &dgO_vec);
  DAVecRestoreArray(BHD->da_newton, u, (void*) &dg_struct);

  /* Compute g's from dg's */
  Zeropad_Function (BHD, dgO, 0.0);
  Zeropad_Function (BHD, dgH, 0.0);
  ComputeH2O_g( gH,  BHD->g_ini[0] , dgH);
  ComputeH2O_g( gO,  BHD->g_ini[1] , dgO);

  //EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH);

  /* Compute right hand side */

  /***********************************************************/
  /* gH */
  /***********************************************************/
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "dgH, ");
  VecCopy(dgH, help);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[0][1], gO, BHD->rhos[1], dgH);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[0][0], gH, BHD->rhos[0], help2);
  VecAXPY(dgH, 1.0, help2);

  VecAXPY(dgH, 1.0, BHD->u2[0][0]);
  /************************************************************/
  /* intra molecular part */
  /************************************************************/
  //goto gH_end;
  Solve_NormalizationH2O_small (BHD, gH, r_HH, gH, tH , help2, ff);
  Compute_dg_H2O_intra_ln(BHD, tH, r_HH, help2);
  VecCopy (help2, ff);          /* FIXME: need that? */
  VecAXPY(dgH, 1.0, help2);
  Solve_NormalizationH2O_small (BHD, gH, r_HO, gO, tO , help2, ff);
  Compute_dg_H2O_intra_ln(BHD, tO, r_HO, help2);
  VecCopy (help2, ff);          /* FIXME: need that? */
  VecAXPY(dgH, 1.0, help2);
  /***********************************************************/
  //gH_end:

/*   ComputeH2O_g( ff,  BHD->g_ini[0] , dgH); */
/*   VecAXPY(gH, -1.0, ff); */
/*   VecCopy(gH, dgH); */

  ImposeLaplaceBoundary (BHD, dgH, BHD->v[0], BHD->v[1]);
  Zeropad_Function (BHD, dgH, 0.0);

  VecAYPX(dgH, -1.0, help);
  //VecPointwiseMult( dgH, dgH, BHD->g_ini[0]);

  /***********************************************************/
  /* gO */
  /***********************************************************/
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "dgO... ");
  VecCopy(dgO, help);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[1][1], gO, BHD->rhos[1], dgO);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[0][1], gH, BHD->rhos[0], help2);
  VecAXPY(dgO, 1.0, help2);

  VecAXPY(dgO, 1.0, BHD->u2[1][1]);
  /************************************************************/
  /* intra molecular part */
  /************************************************************/
  //goto gO_end;
  Solve_NormalizationH2O_small (BHD, gO, r_HO, gH, tH , help2, ff);
  Compute_dg_H2O_intra_ln(BHD, tH, r_HO, help2);
  VecCopy (help2, ff);          /* FIXME: need that? */
  VecAXPY(dgO, 2.0, help2);
  /***********************************************************/
  //gO_end:
/*   ComputeH2O_g( ff,  BHD->g_ini[1] , dgO); */
/*   VecAXPY(gO, -1.0, ff); */
/*   VecCopy(gO, dgO); */

  ImposeLaplaceBoundary (BHD, dgO, BHD->v[0], BHD->v[1]);
  Zeropad_Function (BHD, dgO, 0.0);

  VecAYPX(dgO, -1.0, help);
  //VecPointwiseMult( dgO, dgO, BHD->g_ini[1]);

  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, " done.\n");

  /* Get arrays from PETSC Vectors */
  DAVecGetArray(BHD->da, dgH , &dgH_vec);
  DAVecGetArray(BHD->da, dgO , &dgO_vec);
  DAVecGetArray(BHD->da_newton, f, (void*) &dg_struct);


  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from single Vectors to f */
	  dg_struct[i[2]][i[1]][i[0]].dgH= dgH_vec[i[2]][i[1]][i[0]];
	  dg_struct[i[2]][i[1]][i[0]].dgO= dgO_vec[i[2]][i[1]][i[0]];
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray(BHD->da, dgH , &dgH_vec);
  DAVecRestoreArray(BHD->da, dgO , &dgO_vec);
  DAVecRestoreArray(BHD->da_newton, f, (void*) &dg_struct);

  //VecPointwiseMult(f, f, BHD->pre);
  //VecCopy(f, BHD->pre);

  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "--- Function evaluation finished.\n");

  //WriteH2OSNewtonPlain(BHD, f);

/*   if(counter>250) */
/*     { */
/*       WriteH2OSNewtonPlain(BHD, f);  */
/*       WriteH2OSNewtonSolution(BHD, u); */
/*     }  */
/*   exit(1); */

  return 0;

}

// FIXME: warning of snes cannot be eliminated becaused this would be passed
//      to SNESSetFunction()
static PetscErrorCode ComputeH2OSFunctionFourier(SNES snes, Vec u, Vec f, void *data)
{
  (void) snes;                  /* FIXME: interface obligation? */

  State *BHD;
  H2OSdgF ***dg_struct;
  int i[3], x[3], n[3], index, N3;
  Vec dgH, dgO, gH, gO, help;
  Vec tO, tH, help2, ff;

  static int counter=0;
  counter++;
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "--- Function evaluation starts...\n");

  BHD = (State*) data;
  const ProblemData *PD = BHD->PD;

  gH= BHD->gH;
  gO= BHD->gO;
  dgH= BHD->dgH;
  dgO= BHD->dgO;
  help = BHD->f2;
  help2= BHD->f3;
  tO  =  BHD->f4;
  tH  =  BHD->f4;
  ff = BHD->f;

  Vec uO_fft = BHD->wHO_fft;
  Vec uH_fft = BHD->wHH_fft;

  N3 = PD->N[0]*PD->N[1]*PD->N[2];


  assert (0);                   /* untested */
  /* Get arrays from PETSC Vectors */
  struct {PetscScalar re, im;} ***uH_fft_, ***uO_fft_;
  DAVecGetArray (BHD->da_newtonF, u, (void*) &dg_struct);
  DAVecGetArray (BHD->dc, uH_fft, &uH_fft_);
  DAVecGetArray (BHD->dc, uO_fft, &uO_fft_);

  DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from u to single Vectors */
	  uH_fft_[i[2]][i[1]][i[0]].re = dg_struct[i[2]][i[1]][i[0]].dgHre;
	  uH_fft_[i[2]][i[1]][i[0]].im = dg_struct[i[2]][i[1]][i[0]].dgHim;
	  uO_fft_[i[2]][i[1]][i[0]].re = dg_struct[i[2]][i[1]][i[0]].dgOre;
	  uO_fft_[i[2]][i[1]][i[0]].im = dg_struct[i[2]][i[1]][i[0]].dgOim;
	  index++;
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray (BHD->da_newtonF, u, (void*) &dg_struct);
  DAVecRestoreArray (BHD->dc, uH_fft, &uH_fft_);
  DAVecRestoreArray (BHD->dc, uO_fft, &uO_fft_);

  MatMultTranspose (BHD->fft_mat, uH_fft, dgH);
  MatMultTranspose (BHD->fft_mat, uO_fft, dgO);
  VecScale(dgO, 1.0/N3);
  VecScale(dgH, 1.0/N3);

  Zeropad_Function (BHD, dgO, 0.0);
  Zeropad_Function (BHD, dgH, 0.0);
  ComputeH2O_g( gH,  BHD->g_ini[0] , dgH);
  ComputeH2O_g( gO,  BHD->g_ini[1] , dgO);
  //EnforceNormalizationCondition(BHD, dgO, dgH, gO, gH);

  /* Compute right hand side */

  /***********************************************************/
  /* gH */
  /***********************************************************/
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "dgH, ");
  VecCopy(dgH, help);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[0][1], gO, BHD->rhos[1], dgH);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[0][0], gH, BHD->rhos[0], help2);
  VecAXPY(dgH, 1.0, help2);


  /************************************************************/
  /* intra molecular part */
  /************************************************************/
  //goto gH_end;
  Solve_NormalizationH2O_small (BHD, gH, r_HH, gH, tH , help2, ff);
  Compute_dg_H2O_intra_ln(BHD, tH, r_HH, help2);
  VecCopy (help2, ff);          /* FIXME: need that? */
  VecAXPY(dgH, 1.0, help2);
  Solve_NormalizationH2O_small (BHD, gH, r_HO, gO, tO , help2, ff);
  Compute_dg_H2O_intra_ln(BHD, tO, r_HO, help2);
  VecCopy (help2, ff);          /* FIXME: need that? */
  VecAXPY(dgH, 1.0, help2);
  /***********************************************************/
  //gH_end:

/*   ComputeH2O_g( ff,  BHD->g_ini[0] , dgH); */
/*   VecAXPY(gH, -1.0, ff); */
/*   VecCopy(gH, dgH); */

  ImposeLaplaceBoundary (BHD, dgH, BHD->v[0], BHD->v[1]);
  Zeropad_Function (BHD, dgH, 0.0);

  VecAYPX(dgH, -1.0, help);
  //VecPointwiseMult( dgH, dgH, BHD->g_ini[0]);

  /***********************************************************/
  /* gO */
  /***********************************************************/
  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "dgO... ");
  VecCopy(dgO, help);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[1][1], gO, BHD->rhos[1], dgO);
  Compute_H2O_interS(BHD,
		     BHD->fs_g2_fft[0][1], gH, BHD->rhos[0], help2);
  VecAXPY(dgO, 1.0, help2);
  /************************************************************/
  /* intra molecular part */
  /************************************************************/
  //goto gO_end;
  Solve_NormalizationH2O_small (BHD, gO, r_HO, gH, tH , help2, ff);
  Compute_dg_H2O_intra_ln(BHD, tH, r_HO, help2);
  VecCopy (help2, ff);          /* FIXME: need that? */
  VecAXPY(dgO, 2.0, help2);
  /***********************************************************/
  //gO_end:
/*   ComputeH2O_g( ff,  BHD->g_ini[1] , dgO); */
/*   VecAXPY(gO, -1.0, ff); */
/*   VecCopy(gO, dgO); */

  ImposeLaplaceBoundary (BHD, dgO, BHD->v[0], BHD->v[1]);
  Zeropad_Function (BHD, dgO, 0.0);

  VecAYPX(dgO, -1.0, help);
  //VecPointwiseMult( dgO, dgO, BHD->g_ini[1]);

  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, " done.\n");

  MatMult (BHD->fft_mat, dgH, uH_fft);
  MatMult (BHD->fft_mat, dgO, uO_fft);

  assert (0);                   /* untested */
  /* Get arrays from PETSC Vectors */
  DAVecGetArray(BHD->da_newtonF, f, (void*) &dg_struct);
  DAVecGetArray (BHD->dc, uH_fft, &uH_fft_);
  DAVecGetArray (BHD->dc, uO_fft, &uO_fft_);

  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from single Vectors to f */
	  dg_struct[i[2]][i[1]][i[0]].dgHre= uH_fft_[i[2]][i[1]][i[0]].re;
	  dg_struct[i[2]][i[1]][i[0]].dgHim= uH_fft_[i[2]][i[1]][i[0]].im;
	  dg_struct[i[2]][i[1]][i[0]].dgOre= uO_fft_[i[2]][i[1]][i[0]].re;
	  dg_struct[i[2]][i[1]][i[0]].dgOim= uO_fft_[i[2]][i[1]][i[0]].im;
	  index++;
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray (BHD->da_newtonF, f, (void*) &dg_struct);
  DAVecRestoreArray (BHD->dc, uH_fft, &uH_fft_);
  DAVecRestoreArray (BHD->dc, uO_fft, &uO_fft_);


  if( verbosity>0)
    PetscPrintf(PETSC_COMM_WORLD, "--- Function evaluation finished.\n");

  //WriteH2OSNewtonPlain(BHD, f);

/*   if(counter>250) */
/*     { */
/*       WriteH2OSNewtonPlain(BHD, f);  */
/*       WriteH2OSNewtonSolution(BHD, u); */
/*     }  */
/*   exit(1); */

  return 0;

}
static void WriteH2OSNewtonSolutionF(State *BHD, Vec u)
{
  H2OSdgF ***dg_struct;
  int i[3], x[3], n[3], N3, index;
  Vec dgH, dgO, gH, gO;
  PetscViewer viewer;

  PetscPrintf(PETSC_COMM_WORLD,"Writing files...");

  const ProblemData *PD = BHD->PD;
  gH= BHD->gH;
  gO= BHD->gO;
  dgH= BHD->dgH;
  dgO= BHD->dgO;
  Vec uO_fft = BHD->wHO_fft;
  Vec uH_fft = BHD->wHH_fft;
  N3 = PD->N[0]*PD->N[1]*PD->N[2];

    /* Get arrays from PETSC Vectors */

  struct {PetscScalar re, im;} ***uH_fft_, ***uO_fft_;
  DAVecGetArray (BHD->da_newtonF, u, (void*) &dg_struct);
  DAVecGetArray (BHD->dc, uH_fft, &uH_fft_);
  DAVecGetArray (BHD->dc, uO_fft, &uO_fft_);


  DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));
  index=0;
  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from u to single Vectors */
	  uH_fft_[i[2]][i[1]][i[0]].re = dg_struct[i[2]][i[1]][i[0]].dgHre;
	  uH_fft_[i[2]][i[1]][i[0]].im = dg_struct[i[2]][i[1]][i[0]].dgHim;
	  uO_fft_[i[2]][i[1]][i[0]].re = dg_struct[i[2]][i[1]][i[0]].dgOre;
	  uO_fft_[i[2]][i[1]][i[0]].im = dg_struct[i[2]][i[1]][i[0]].dgOim;
	  index++;
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray (BHD->da_newtonF, u, (void*) &dg_struct);
  DAVecRestoreArray (BHD->dc, uH_fft, &uH_fft_);
  DAVecRestoreArray (BHD->dc, uO_fft, &uO_fft_);

  MatMultTranspose (BHD->fft_mat, uH_fft, dgH);
  MatMultTranspose (BHD->fft_mat, uO_fft, dgO);
  VecScale(dgO, 1.0/N3);
  VecScale(dgH, 1.0/N3);


  /* Copmute g's from dg's */
  ComputeH2O_g( gH,  BHD->g_ini[0] , dgH);
  ComputeH2O_g( gO,  BHD->g_ini[1] , dgO);

  /*************************************/
  /* output */
  /* g_H */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecH.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecH.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gH,viewer);
  PetscViewerDestroy(viewer);
  /* g_O */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"vecO.m",&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecO.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gO,viewer);
  PetscViewerDestroy(viewer);
  PetscPrintf(PETSC_COMM_WORLD,"done\n");
  /************************************/


}


static void WriteH2OSNewtonSolution(State *BHD, Vec u)
{
  H2OSdg ***dg_struct;
  PetscScalar ***dgH_vec, ***dgO_vec;
  int i[3], x[3], n[3];
  Vec dgH, dgO, gH, gO;
  PetscViewer viewer;
  static int count=0;
  char nameO[20], nameH[20];

  PetscPrintf(PETSC_COMM_WORLD,"Writing files...");

  //count++;
  sprintf(nameO, "vecO%d.m", count-1);
  sprintf(nameH, "vecH%d.m", count-1);


  gH= BHD->gH;
  gO= BHD->gO;
  dgH= BHD->dgH;
  dgO= BHD->dgO;


  /* Get arrays from PETSC Vectors */
  DAVecGetArray(BHD->da, dgH, &dgH_vec);
  DAVecGetArray(BHD->da, dgO, &dgO_vec);
  DAVecGetArray(BHD->da_newton, u, (void*) &dg_struct);

  DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  /* Copy from u to single Vectors */
	  dgH_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgH;
	  dgO_vec[i[2]][i[1]][i[0]] = dg_struct[i[2]][i[1]][i[0]].dgO;
	}
  /* Restore arrays from PETSC Vectors */
  DAVecRestoreArray(BHD->da, dgH, &dgH_vec);
  DAVecRestoreArray(BHD->da, dgO, &dgO_vec);
  DAVecRestoreArray(BHD->da_newton, u, (void*) &dg_struct);

  /* Copmute g's from dg's */
  ComputeH2O_g( gH,  BHD->g_ini[0] , dgH);
  ComputeH2O_g( gO,  BHD->g_ini[1] , dgO);

  /*************************************/
  /* output */
  /* g_H */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameH,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecH.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gH,viewer);
  PetscViewerDestroy(viewer);
  /* g_O */
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,nameO,&viewer);
  //PetscViewerBinaryOpen(PETSC_COMM_WORLD,"vecO.m",FILE_MODE_WRITE,&viewer);
  PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
  VecView(gO,viewer);
  PetscViewerDestroy(viewer);
  PetscPrintf(PETSC_COMM_WORLD,"done\n");
  /************************************/


}







/* apply preconditioner matrix: diagonal scaling */
static PetscErrorCode ComputePreconditioner_H2OS(void *data, Vec x, Vec y)
{
  State *BHD;
  PetscErrorCode ierr;

  BHD = (State*) data;
  //ierr=VecPointwiseMult(y,BHD->pre,x);

  VecAbs(BHD->pre);
  VecShift(BHD->pre, 1.0);
  ierr = VecPointwiseDivide(y, x, BHD->pre);

/*   VecView(x,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);  */

  return ierr;
}

#include "petscmat.h"

Vec BGY3d_SolveNewton_H2OS(const ProblemData *PD, Vec g_ini)
{
  SNES snes;
  KSP ksp;
  PC  pc;
  State *BHD;
  Vec u, f, b, v1, v2;
  real damp;
  // Mat M, A;
  PetscTruth flg;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Newton ...\n");


  BHD = BGY3dH2OData_Newton_malloc(PD);

  /* Get damp_start from command line*/
  const real damp_start = PD->damp;

  DACreateGlobalVector(BHD->da_newton, &f);
  DACreateGlobalVector(BHD->da_newton, &u);
  DACreateGlobalVector(BHD->da_newton, &b);
  DACreateGlobalVector(BHD->da_newton, &v1);
  DACreateGlobalVector(BHD->da_newton, &v2);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD->da, BHD->PD, &BHD->M, &BHD->ksp);
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
    PCShellSetApply(pc,ComputePreconditioner_H2OS);
    PCShellSetContext(pc,BHD);
  } else
    /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
    PCSetType( pc, PCJACOBI);
  /* set function */
  SNESSetFunction(snes, f, ComputeH2OSFunction, (void*)BHD);

  /* runtime options will override default parameters */
  SNESSetFromOptions(snes);

  /* set initial guess */
  VecSet(u, 0.0);
  //VecSetRandom_H2O(u, 0.5);

/*   MatCreateSNESMF(snes, u, &M); */
/*   MatCreateSeqDense(PETSC_COMM_WORLD,3*PD->N3,3*PD->N3,PETSC_NULL,&A); */
/*   SNESDefaultComputeJacobian(snes,u,&A,&A,SAME_NONZERO_PATTERN,BHD); */

  for(damp=damp_start; damp<=1; damp+= 0.01)
    {

      RecomputeInitialFFTs(BHD, (damp), 1.0);
      RecomputeInitialSoluteData(BHD, (damp), 1.0);
      ImposeLaplaceBoundary (BHD, BHD->g_ini[0], BHD->v[0], BHD->v[1]);
      ImposeLaplaceBoundary (BHD, BHD->g_ini[1], BHD->v[0], BHD->v[1]);
      Zeropad_Function (BHD, BHD->g_ini[0], 0.0);
      Zeropad_Function (BHD, BHD->g_ini[1], 0.0);
      InitializePreconditioner( BHD);

      /* solve problem */
      SNESSolve(snes, PETSC_NULL, u);

      /* Get Solution */
      //SNESGetSolution(snes, &u);
      SNESComputeFunction(snes,u,f);

      /* Write out solution */
      WriteH2OSNewtonSolution(BHD, u);
      //WriteH2OSNewtonPlain(BHD, f);
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



Vec BGY3d_SolveNewton_H2OSF(const ProblemData *PD, Vec g_ini)
{
  SNES snes;
  KSP ksp;
  PC  pc;
  State *BHD;
  Vec u, f, b;
  real damp;
  // Mat M, A;
  PetscTruth flg;

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Newton ...\n");


  BHD = BGY3dH2OData_Newton_malloc(PD);

  /* Get damp_start from command line*/
  const real damp_start = PD->damp;

  DACreateGlobalVector(BHD->da_newtonF, &f);
  DACreateGlobalVector(BHD->da_newtonF, &u);
  DACreateGlobalVector(BHD->da_newtonF, &b);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD->da, BHD->PD, &BHD->M, &BHD->ksp);
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
    PCShellSetApply(pc,ComputePreconditioner_H2OS);
    PCShellSetContext(pc,BHD);
  } else
    /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
    PCSetType( pc, PCJACOBI);
  /* set function */
  SNESSetFunction(snes, f, ComputeH2OSFunctionFourier, (void*)BHD);

  /* runtime options will override default parameters */
  SNESSetFromOptions(snes);

  /* set initial guess */
  VecSet(u, 0.0);
  //VecSetRandom_H2O(u, 0.5);

/*   MatCreateSNESMF(snes, u, &M); */
/*   MatCreateSeqDense(PETSC_COMM_WORLD,3*PD->N3,3*PD->N3,PETSC_NULL,&A); */
/*   SNESDefaultComputeJacobian(snes,u,&A,&A,SAME_NONZERO_PATTERN,BHD); */

  for(damp=damp_start; damp<=1; damp+= 0.1)
    {

      RecomputeInitialFFTs(BHD, SQR(damp), 1.0);
      RecomputeInitialSoluteData(BHD, SQR(damp), 1.0);
      ImposeLaplaceBoundary (BHD, BHD->g_ini[0], BHD->v[0], BHD->v[1]);
      ImposeLaplaceBoundary (BHD, BHD->g_ini[1], BHD->v[0], BHD->v[1]);
      Zeropad_Function(BHD, BHD->g_ini[0], 0.0);
      Zeropad_Function(BHD, BHD->g_ini[1], 0.0);
      InitializePreconditioner( BHD);

      /* solve problem */
      SNESSolve(snes, PETSC_NULL, u);

      /* Get Solution */
      //SNESGetSolution(snes, &u);
      SNESComputeFunction(snes,u,f);

      /* Write out solution */
      WriteH2OSNewtonSolutionF(BHD, u);
      //WriteH2OSNewtonPlain(BHD, f);
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

static void RecomputeInitialSoluteData(State *BHD, real damp, real damp_LJ)
{
  DA da;
  PetscScalar ***gHini_vec, ***gOini_vec;
  // PetscScalar ***ucH_vec, ***ucO_vec;
  PetscScalar ***(fHl_vec[3]),***(fOl_vec[3]);
  real r[3], r_s, h[3], interval[2], beta;
  int x[3], n[3], i[3];
  real epsilonH, epsilonO, epsilonHO;
  real sigmaH, sigmaO, sigmaHO;
  real q2H, q2O, q2HO;

  assert (0);                   /* See FIXME below! */

  epsilonH = BHD->LJ_paramsH[0];
  epsilonO = BHD->LJ_paramsO[0]; /* FIXME: how comes it is unused? */
  epsilonHO = BHD->LJ_paramsHO[0];

  sigmaH = BHD->LJ_paramsH[1];
  sigmaO = BHD->LJ_paramsO[1];  /* FIXME: how comes it is unused? */
  sigmaHO = BHD->LJ_paramsHO[1];

  q2H = BHD->LJ_paramsH[2];
  q2O = BHD->LJ_paramsO[2];     /* FIXME: how comes it is unused? */
  q2HO = BHD->LJ_paramsHO[2];


  const ProblemData *PD = BHD->PD;
  da = BHD->da;

  PetscPrintf(PETSC_COMM_WORLD,"Recomputing solute data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);


  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  beta = PD->beta;

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);



  VecSet(BHD->g_ini[0], 0.0);
  VecSet(BHD->g_ini[1], 0.0);
  VecSet(BHD->gHO_ini, 0.0);
  VecSet(BHD->u2[0][0], 0.0);
  VecSet(BHD->u2[1][1], 0.0);
  VecSet(BHD->u2[0][1], 0.0);
  FOR_DIM
    {
      VecSet(BHD->F_l[0][0][dim],0.0);
      VecSet(BHD->F_l[1][1][dim],0.0);
      VecSet(BHD->F_l[0][1][dim],0.0);
    }


  DAVecGetArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecGetArray(da, BHD->g_ini[1], &gOini_vec);

/*   DAVecGetArray(da, BHD->u2[0][0], &ucH_vec); */
/*   DAVecGetArray(da, BHD->u2[1][1], &ucO_vec); */
  FOR_DIM
    {
      DAVecGetArray(da, BHD->F_l[0][0][dim], &fHl_vec[dim]);
      DAVecGetArray(da, BHD->F_l[1][1][dim], &fOl_vec[dim]);
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
/*        gHini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsH); */
/*        gOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */
/*        gHOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */

           /* Coulomb */
/*        gHini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsH); */
/*        gOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsO); */
/*        gHOini_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb( r_s, BHD->LJ_paramsHO); */

/*        ucH_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsHO); */
/*        ucO_vec[i[2]][i[1]][i[0]] +=  */
/*          damp * beta* Coulomb_long( r_s, BHD->LJ_paramsO); */



        }

  DAVecRestoreArray(da, BHD->g_ini[0], &gHini_vec);
  DAVecRestoreArray(da, BHD->g_ini[1], &gOini_vec);
/*   DAVecRestoreArray(da, BHD->u2[0][0], &ucH_vec); */
/*   DAVecRestoreArray(da, BHD->u2[1][1], &ucO_vec); */

  VecCopy( BHD->u2[0][1], BHD->u2[1][1]);
}
