/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-vec.h"          /* bgy3d_vec_duplicate() */
#include "bgy3d-snes.h"         /* Function, Solver */

/*
  For solving HNC equation with Newton. Except of x everything else in
  the  closure context is  considered read  only input  or intemediate
  terms depending  on x.  That  is when looking for  total correlation
  function h,  the direct correlation  function should be fixed  (or a
  function of h).
*/
void bgy3d_snes_newton (const ProblemData *PD, void *ctx, Function F, Vec x)
{
  (void) PD;             /* FIXME: convergence criteria are ignored */

  /* Create the snes environment */
  SNES snes;
  SNESCreate (PETSC_COMM_WORLD, &snes);

  KSP ksp;
  SNESGetKSP (snes, &ksp);

  PC pc;
  KSPGetPC (ksp, &pc);

  /* set rtol, atol, dtol, maxits */
  KSPSetTolerances (ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);

  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType (snes, SNESLS);

  /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
  PCSetType (pc, PCNONE);

  /* SNES needs a place to store residual: */
  Vec r = bgy3d_vec_duplicate (x);

  /* SNES functions should obey this interface: */
  PetscErrorCode F1 (SNES snes, Vec x, Vec r, void *ctx)
  {
    (void) snes;                /* unused */
    F (ctx, x, r);              /* assumes ctx is a Context* */
    return 0;
  }
  SNESSetFunction (snes, r, F1, ctx); /* Pass Context* as ctx */

  /* set atol, rtol, stol , its, fct. eval. */
  // SNESSetTolerances (snes, 5.0e-2, 1.0e-5, 1.0e-4 , 50, 10000);
  // SNESSetTolerances (snes, 5.0e-2, 1.0e-5, PD->norm_tol, 50, 10000);

  /*
    Runtime  options will  override default  parameters.   FIXME: note
    that the  call to SNESSetJacobian()  is missing here.   It appears
    that  one has to  request a  "matrix-free" approximation  from the
    command line  with "-snes_mf". Otherwise the  next call terminates
    with an error message saying "Matrix must be set first"!
  */
  SNESSetFromOptions (snes);

  /* Solve  problem F(x)  = 0.  PETSC_NULL indicates  that the  rhs is
     0: */
  SNESSolve (snes, PETSC_NULL, x);

  /* Write  out  solution.   FIXME:   and  what  was  the  purpose  of
     SNESSolve()? */
  // SNESGetSolution (snes, &x);

  bgy3d_vec_destroy (&r);

  SNESDestroy (snes);
}

void bgy3d_snes_picard (const ProblemData *PD, void *ctx, Function F, Vec x)
{
  /* Mixing parameter */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* Convergence threshold: */
  const real norm_tol = PD->norm_tol;

  /* A place to store residual: */
  Vec dx = bgy3d_vec_duplicate (x);

  /* Find an x such that dx as returned by F (ctx, x, dx) is zero: */
  for (int k = 0; k < max_iter; k++)
    {
      F (ctx, x, dx);

      /* Simple mixing: x = lambda * x + (1 - lambda) * x_old */
      VecAXPY (x, lambda, dx);

      const real norm = bgy3d_vec_norm (dx);

      PetscPrintf (PETSC_COMM_WORLD, "%03d: norm of difference: %e\t%f\n",
                   k + 1, norm, lambda);

      if (norm < norm_tol)
        break;
    }
  bgy3d_vec_destroy (&dx);
}
