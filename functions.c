/*==========================================================*/
/*  $Id: functions.c,v 1.8 2006-10-05 16:37:44 jager Exp $ */
/*==========================================================*/


#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "functions.h"

static PetscErrorCode ComputeVec_F(SNES snes, Vec g, Vec f, void *pa);

typedef struct BGY3dVecStruct
{
  BGY3dParameter params[3];
  Vec fl[3];

} *BGY3dParameterVec;

static BGY3dParameterVec BGY3dParameterVec_malloc(ProblemData *PD);
static void BGY3dParameterVec_free(BGY3dParameterVec par_vec);
static void CreateInitialGuess_vec(BGY3dParameterVec par_vec, Vec g);

Vec BGY3d_vec_solve(ProblemData *PD, Vec g_ini, int vdim)
{
  Vec F, g;
  BGY3dParameterVec par_vec;
  SNES snes;
  KSP ksp;
  PC  pc;
  PetscTruth flg;

  par_vec = BGY3dParameterVec_malloc(PD);

  /* Create global vectors */
  VecDuplicate(par_vec->params[0]->x, &F);
  VecDuplicate(par_vec->params[0]->x, &g);

  CreateInitialGuess_vec(par_vec, g);

  /* Create the snes environment */
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp, &pc);

  /* set rtol, atol, dtol, maxits */
  // KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  KSPSetTolerances(ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);
  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType(snes, SNESLS);

  flg = bgy3d_getopt_test ("-user_precond");
  if (flg) { /* user-defined precond */
    /* Set user defined preconditioner */
    PCSetType(pc,PCSHELL);
    PCShellSetApply(pc,Compute_Preconditioner);
    PCShellSetContext(pc,par_vec->params[0]);
  } else
    /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
    PCSetType( pc, PCNONE);

  /* ComputeVec_F(snes, g, F, (void*) par_vec); */
/*   Global2Local(par_vec, gl, F); */
/*   VecView(gl[0],PETSC_VIEWER_STDERR_WORLD); */
/*   exit(1); */

  SNESSetFunction(snes, F, ComputeVec_F, par_vec);

  /* runtime options will override default parameters */
  SNESSetFromOptions(snes);

  /* solve problem */
  SNESSolve(snes, PETSC_NULL, g);



  /* write out solution */
  SNESGetSolution(snes, &g);



  /* free stuff */
  BGY3dParameterVec_free(par_vec);

  VecDestroy(F);
  SNESDestroy(snes);

  return g;

}

static BGY3dParameterVec BGY3dParameterVec_malloc(ProblemData *PD)
{
  BGY3dParameterVec par_vec;

  par_vec = (BGY3dParameterVec) malloc(sizeof(*par_vec));



  /* Create parameter structs */
  FOR_DIM
    par_vec->params[dim] = BGY3dParameter_malloc(PD, dim);

  FOR_DIM
    VecDuplicate(par_vec->params[dim]->x, &(par_vec->fl[dim]));

  return par_vec;
}


static void BGY3dParameterVec_free(BGY3dParameterVec par_vec)
{
  FOR_DIM
    {
      VecDestroy(par_vec->fl[dim]);
      BGY3dParameter_free(par_vec->params[dim]);
    }

  free(par_vec);
}


static void CreateInitialGuess_vec(BGY3dParameterVec par_vec, Vec g)
{

  CreateInitialGuess(par_vec->params[0], g);
  //VecSet(g,0.0);
}





static PetscErrorCode ComputeVec_F(SNES snes, Vec g, Vec f, void *pa)
{
  BGY3dParameterVec par_vec;

  par_vec = (BGY3dParameterVec) pa;

  for(int dim = 0; dim < 3; dim++)
    Compute_F(snes, g, par_vec->fl[dim],
		       (void *) par_vec->params[dim]);

  VecSet(f,0.0);
  FOR_DIM
    {
      //VecAbs(par_vec->fl[dim]);
      VecPointwiseMult(par_vec->fl[dim],par_vec->fl[dim],par_vec->fl[dim]);
      VecAXPY(f,1.0,par_vec->fl[dim]);
    }

  return 0;
}
