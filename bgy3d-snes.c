/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-vec.h"          /* bgy3d_vec_duplicate() */
#include "bgy3d-getopt.h"       /* bgy3d_getopt_string() */
#include "bgy3d-snes.h"         /* Function, Solver */


void bgy3d_snes_default (const ProblemData *PD, void *ctx, Function F, Vec x)
{
  char solver[20] = "newton";
  bgy3d_getopt_string ("--snes-solver", solver, sizeof solver);

  if (strcmp (solver, "newton") == 0)
    bgy3d_snes_newton (PD, ctx, F, x);
  else if (strcmp (solver, "picard") == 0)
    bgy3d_snes_picard (PD, ctx, F, x);
  else if (strcmp (solver, "jager") == 0)
    bgy3d_snes_jager (PD, ctx, F, x);
  else
    {
      PetscPrintf (PETSC_COMM_WORLD, "No such SNES solver: %s\n", solver);
      exit (1);
    }
}


/*
  For solving HNC equation with Newton. Except of x everything else in
  the  closure context is  considered read  only input  or intemediate
  terms depending  on x.  That  is when looking for  total correlation
  function h,  the direct correlation  function should be fixed  (or a
  function of h).
*/
void bgy3d_snes_newton (const ProblemData *PD, void *ctx, Function F, Vec x)
{
  /* Create the snes environment */
  SNES snes;
  SNESCreate (PETSC_COMM_WORLD, &snes);

  /* SNES needs a place to store residual: */
  Vec r = bgy3d_vec_duplicate (x);

  /* SNES form-functions should obey this interface: */
  PetscErrorCode F1 (SNES snes, Vec x, Vec r, void *ctx)
  {
    (void) snes;                /* unused */
    F (ctx, x, r);              /* assumes ctx is a Context* */
    return 0;
  }
  SNESSetFunction (snes, r, F1, ctx); /* Pass Context* as ctx */

  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType (snes, SNESLS);

  {
    /*
      This has the same effect as specifying "-snes_mf" in the command
      line.  We  do it  explicitly here in  order not to  require that
      switch from the user. MF stays for Matrix-Free.

      The  idea  is  that  matrix-vector products  with  Jacobian  are
      computed   by   numerical   differentiation  of   the   original
      form-function.  The  form-function has to  already be associated
      with the SNES object, see above.
    */
    Mat J;                      /* I guess we need to destroy it? */
    MatCreateSNESMF (snes, &J);

    /*
      Petsc provides a placeholder for this case:

        PetscErrorCode
        MatMFFDComputeJacobian (SNES, Vec, Mat*, Mat*, MatStructure*, void*)

      The last argument is a user context for jacobian evaluation:
    */
    SNESSetJacobian (snes, J, J, MatMFFDComputeJacobian, NULL);
    MatDestroy (J);         /* I hope SNES saved a ref to that Mat? */
  }

  /* set atol, rtol, stol , its, fct. eval. */
  {
    real atol, rtol, stol;
    int max_it, max_funcs;

    /* Only max-iter  changed. Sometimes it  is convenient to  let the
       solver do just one, or just a few iterations: */
    SNESGetTolerances (snes, &atol, &rtol, &stol, &max_it, &max_funcs);
    SNESSetTolerances (snes, atol, rtol, stol, PD->max_iter, max_funcs);
    // SNESSetTolerances (snes, 5.0e-2, 1.0e-5, 1.0e-4 , 50, 10000);
  }

  {
    /* This linear equation solver is  most probably used to solve the
       linear equations with Jacobian matrix: */
    KSP ksp;
    SNESGetKSP (snes, &ksp);    /* no need to destroy, apparently */

    /* set rtol, atol, dtol, maxits */
    KSPSetTolerances (ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);

    /*
      Though the manual says:

        "The  matrix-free  variant is  allowed  only  when the  linear
         systems are solved by an iterative method in combination with
         no preconditioning (PCNONE or -pc_type none), a user-provided
         preconditioner  matrix,  or  a  user-provided  preconditioner
         shell (PCSHELL, discussed in Section 4.4); that is, obviously
         matrix-free methods cannot  be used if a direct  solver is to
         be employed."

      it seems that commenting the following lines does not affect the
      behaviour.  FIXME:  Didnt we  set the preconditioner  with third
      argument of SNESSetJacobian()?
    */
    PC pc;
    KSPGetPC (ksp, &pc);

    /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
    PCSetType (pc, PCNONE);
  }

  /*
    Runtime options  will override  default parameters.  Note  that if
    SNESSetJacobian() is not called,  one has to request a matrix-free
    approximation  from the command  line with  "-snes_mf".  Otherwise
    the next call terminates with an error message saying "Matrix must
    be set first"!
  */
  SNESSetFromOptions (snes);

  /* Solve  problem F(x)  = 0.  PETSC_NULL indicates  that the  rhs is
     0: */
  SNESSolve (snes, PETSC_NULL, x);

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

void bgy3d_snes_jager (const ProblemData *PD, void *ctx, Function F, Vec x)
{
  /* Mixing parameter */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* Convergence threshold: */
  const real norm_tol = PD->norm_tol;

  /* A place to store residual: */
  Vec dx = bgy3d_vec_duplicate (x);

  /* Not sure if 0.0 as inital value is right. */
  real norm_old = 0.0;

  /* Find an x such that dx as returned by F (ctx, x, dx) is zero: */
  const real a0 = lambda;
  real a1 = lambda;             /* loop-local variable */
  for (int iter = 0, mycount = 0, upwards = 0; iter < max_iter; iter++)
    {
      /* Calculate residual: */
      F (ctx, x, dx);

      const real norm = bgy3d_vec_norm (dx);

      /*
        Most  of  the  logic  below  is  to  control  how  the  mixing
        coefficient  "a" changes  from interation  to  iterations also
        dependign on the behaviour of the residual norm.
      */
      const int nth = 10;
      /*
        "a  = a1"  is taken  in  iteration 0,  10, 20,  etc.  "a1"  is
        modified during the loop.

        "a = a0" is taken in iterations 1-9, 11-19, etc.  "a0" remains
        unchanged during the loop.

        Note that in the first iteration a1 == a0.
      */

      /* Every nth  iteration, raise the mixing  coefficients just one
         time: */
      const real a = (iter % nth == 0) ? a1 : a0;

      /* Simple mixing: x = a * x + (1 - a) * x_old */
      VecAXPY (x, a, dx);

      /* Fancy step  size control. FIXME:  weired logic. Code  used to
         check if the norm went up: */
      const bool up = norm > norm_old;

      /* That was the only place comparing to norm_old: */
      norm_old = norm;

      mycount++;

      if (iter % nth != 1 && up) /* not in the nth + 1 iteration */
        upwards = 1;
      else if (iter > 2 * nth && iter % nth == 1 && upwards == 0 && up)
        {
          /* In the  nth + 1 iteration,  if the norm  went up decrease
             the mixing: */
          a1 = MAX (a1 / 2.0, a0);
          mycount = 0;
        }
      else
        upwards = 0;

      /* Scale the coefficient  "a1" up by a factor,  but make sure it
         is not above 1.0. Reset mycount. */
      if (mycount > 2 * nth)
        {
          a1 = MIN (a1 * 2.0, 1.0);
          mycount = 0;
        }
      /* otherwise leave "a1" and "mycount" unchanged */

      PetscPrintf (PETSC_COMM_WORLD, "%03d: norm of difference: %e\t%f",
                   iter + 1, norm, a);
      PetscPrintf (PETSC_COMM_WORLD, " count=%3d upwards=%1d", mycount, upwards);
      PetscPrintf (PETSC_COMM_WORLD, "\n");

      /* Exit when residual norm does not exceed norm_tol: */
      if (norm <= norm_tol)
        {
          PetscPrintf (PETSC_COMM_WORLD,
                       "norm %e <= %e (norm-tol) in iteration %d < %d (max-iter)\n",
                       norm, norm_tol, iter + 1, max_iter);
          break;
        }
    } /* for (iter = ... ) */
  bgy3d_vec_destroy (&dx);
}
