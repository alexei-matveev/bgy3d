/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-vec.h"          /* vec_duplicate() */
#include "bgy3d-getopt.h"       /* bgy3d_getopt_string() */
#include "bgy3d-snes.h"         /* VectorFunc, ArrayFunc */


void bgy3d_snes_default (const ProblemData *PD, void *ctx, VectorFunc F, Vec x)
{
  char solver[20] = "newton";
  bgy3d_getopt_string ("--snes-solver", sizeof solver, solver);

  if (strcmp (solver, "newton") == 0)
    bgy3d_snes_newton (PD, ctx, F, x);
  else if (strcmp (solver, "picard") == 0)
    bgy3d_snes_picard (PD, ctx, F, x);
  else if (strcmp (solver, "jager") == 0)
    bgy3d_snes_jager (PD, ctx, F, x);
  else if (strcmp (solver, "trial") == 0)
    bgy3d_snes_trial (PD, ctx, F, x);
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
  function h, the direct correlation function should be fixed (or be a
  function of h).
*/
void bgy3d_snes_newton (const ProblemData *PD, void *ctx, VectorFunc F, Vec x)
{
  /* Create the snes environment */
  SNES snes;
  SNESCreate (PETSC_COMM_WORLD, &snes);

  /* SNES needs a place to store residual: */
  local Vec r = vec_duplicate (x);

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
    MatDestroy (&J);         /* I hope SNES saved a ref to that Mat? */
  }

  /* set atol, rtol, stol , its, fct. eval. */
  {
    real atol, rtol, stol;
    int max_it, max_funcs;

    /* Only max-iter  changed. Sometimes it  is convenient to  let the
       solver do just one, or just a few iterations: */
    SNESGetTolerances (snes, &atol, &rtol, &stol, &max_it, &max_funcs);
    SNESSetTolerances (snes, atol, PD->norm_tol, stol, PD->max_iter, max_funcs);
    // SNESSetTolerances (snes, 5.0e-2, 1.0e-5, 1.0e-4 , 50, 10000);
  }

  {
    /* This linear equation solver is  most probably used to solve the
       linear equations with Jacobian matrix: */
    KSP ksp;
    SNESGetKSP (snes, &ksp);    /* no need to destroy, apparently */

    /* set rtol, atol, dtol, maxits */
    {
      real rtol, abstol, dtol;
      int maxits;
      KSPGetTolerances (ksp, &rtol, &abstol, &dtol, &maxits);
      /*
        Defaults are at the mercy  of the library: rtol = 1e-5, abstol
        = 1e-50, dtol = 1e+4, maxits = 10000.

        Note that  these convergence criteria control  the accuracy of
        the Jacobian sub-problems used to determine the "direction" of
        the  line search  for the  non-linear equation  system.  Thus,
        these convergence  criteria will  probably affect the  rate of
        SNES convergence, but not  the accuracy of SNES solution. With
        the default  convergence criteria the impression  was that the
        solver  wastes  unnecessarily much  time  refining  δx in  the
        equivalent of linear equation J δx = -r in the region far from
        SNES convergence.  Therefore,  the somewhat drastic factor 100
        here.   With the defaults  above the  code below  requires the
        norm of  the residual of the *linear*  Jacobian sub-problem to
        be reduced  by at  least factor  1000 in order  to be  used as
        search direction:
      */
      KSPSetTolerances (ksp, 100 * rtol, abstol, 10 * dtol, maxits / 10);
    }

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

  /*
    It looks like SNESGetSolution() is only of any value for callbacks
    that need  to extract intermediate solution from  the SNES object.
    Here it is  fully redundant.  Do not vec_destroy  (&y), check this
    assert out. Appears to hold  even when SNES does not converge, say
    due to iteration limit:
  */
  {
    Vec y;
    SNESGetSolution (snes, &y); /* borrow Vec */
    assert (x == y);
  }

  vec_destroy (&r);

  SNESDestroy (&snes);
}

void bgy3d_snes_picard (const ProblemData *PD, void *ctx, VectorFunc F, Vec x)
{
  /* Mixing parameter */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* Convergence threshold: */
  const real norm_tol = PD->norm_tol;

  /* A place to store residual: */
  local Vec dx = vec_duplicate (x);

  /* Find an x such that dx as returned by F (ctx, x, dx) is zero: */
  for (int k = 0; k < max_iter; k++)
    {
      F (ctx, x, dx);

      /* Simple mixing: x = lambda * x + (1 - lambda) * x_old */
      VecAXPY (x, lambda, dx);

      const real norm = vec_norm (dx);

      if (verbosity > 0)
        PetscPrintf (PETSC_COMM_WORLD, " # %03d: norm of difference: %e\t%f\n",
                     k + 1, norm, lambda);

      if (norm < norm_tol)
        break;
    }
  vec_destroy (&dx);
}


void bgy3d_snes_trial (const ProblemData *PD, void *ctx, VectorFunc F, Vec x)
{
  /* First do a few "slow" iterations: */
  {
    ProblemData pd = *PD;       /* modify a copy, not the original */
    pd.max_iter = 100;
    pd.lambda = 0.02;
    bgy3d_snes_picard (&pd, ctx, F, x);
  }
  /* Then continue with Newton: */
  bgy3d_snes_newton (PD, ctx, F, x);
}


void bgy3d_snes_jager (const ProblemData *PD, void *ctx, VectorFunc F, Vec x)
{
  /* Mixing parameter */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* Convergence threshold: */
  const real norm_tol = PD->norm_tol;

  /* A place to store residual: */
  local Vec dx = vec_duplicate (x);

  /* Not sure if 0.0 as inital value is right. */
  real norm_old = 0.0;

  /* Find an x such that dx as returned by F (ctx, x, dx) is zero: */
  const real a0 = lambda;
  real a1 = lambda;             /* loop-local variable */
  for (int iter = 0, mycount = 0, upwards = 0; iter < max_iter; iter++)
    {
      /* Calculate residual: */
      F (ctx, x, dx);

      const real norm = vec_norm (dx);

      /*
        Most  of  the  logic  below  is  to  control  how  the  mixing
        coefficient  "a" changes  from interation  to  iterations also
        dependign on the behaviour of the residual norm.
      */
      const int nth = 10;
      /*
        "a  = a1"  is taken  in  iteration 0,  10, 20,  etc.  "a1"  is
        modified during the loop (moves).

        "a = a0" is taken in iterations 1-9, 11-19, etc.  "a0" remains
        unchanged during the loop (annealing).

        Note that in the first iteration a1 == a0.
      */

      /* Every  nth  iteration make  a  "move"  ---  raise the  mixing
         coefficients just one time: */
      const real a = (iter % nth == 0) ? a1 : a0;

      /* Simple mixing: x = a * x + (1 - a) * x_old */
      VecAXPY (x, a, dx);

      /* Fancy step  size control. FIXME:  weired logic. Code  used to
         check if the norm went up: */
      const bool up = norm > norm_old;

      /* That was the only place comparing to norm_old: */
      norm_old = norm;

      /* Measure time since the last change of "a1": */
      mycount++;

      /* 1) Watching "annealing": */
      if (iter % nth != 1)      /* not in the nth + 1 iteration ... */
        {
          if (up)               /* if the norm went up ... */
            upwards = 1;        /* we are not even near convergence. */
          else
            upwards = 0;        /* otherwise we might be! */
        }

      /* 2) Watching "moves": */
      if (iter % nth == 1)      /* in the  nth + 1 iteration ... */
        if (up)                 /* if the norm went up ... */
          if (upwards == 0)     /* but otherwise it did not ... */
            if (iter > 2 * nth) /* and we did not just start ... */
              {
                /*
                  The "move"  was too dangerous.   Decrease the mixing
                  for the "moves", but not below that of "annealing":
                */
                a1 = MAX (a1 / 2.0, a0);
                mycount = 0;
              }

      /*
        3) Make more  couragous "moves"  as  time  passes.   Scale the
           coefficient "a1"  up by a factor,  but make sure  it is not
           above 1.0. Reset mycount.
      */
      if (mycount > 2 * nth)
        {
          a1 = MIN (a1 * 2.0, 1.0);
          mycount = 0;
        }
      /* otherwise leave "a1" and "mycount" unchanged */

      if (verbosity > 0)
        {
          PetscPrintf (PETSC_COMM_WORLD, " # %03d: norm of difference: %e\t%f",
                       iter + 1, norm, a);
          PetscPrintf (PETSC_COMM_WORLD, " count=%3d upwards=%1d", mycount, upwards);
          PetscPrintf (PETSC_COMM_WORLD, "\n");
        }

      /* Exit when residual norm does not exceed norm_tol: */
      if (norm <= norm_tol)
        {
          if (verbosity > 0)
            PetscPrintf (PETSC_COMM_WORLD,
                         " # norm %e <= %e (norm-tol) in iteration %d < %d (max-iter)\n",
                         norm, norm_tol, iter + 1, max_iter);
          break;
        }
    } /* for (iter = ... ) */
  vec_destroy (&dx);
}


/*
  A solver  for an untyped (array) form-function  ArrayFunc.  Finds an
  x_[n] such that dx_[n] as returned  by ArrayFunc f (ctx, n, x_, dx_)
  is zero.  May be eventually used  from Fortran, so  if ever changing
  the interface  update the  corresponding interface block  in Fortran
  sources.
*/
void rism_snes (void *ctx, ArrayFunc f, int n, real x_[n])
{
  local Vec x = vec_from_array (n, x_);

  /* Implements VectorFunc interface: */
  void F (void *ctx, Vec y, Vec dy)
  {
    local real *y_ = vec_get_array (y);
    local real *dy_ = vec_get_array (dy);

    /* Implements ArrayFunc interface: */
    f (ctx, n, y_, dy_);

    vec_restore_array (y, &y_);
    vec_restore_array (dy, &dy_);
  }

  /* The result depends on the command line and affects the solver: */
  ProblemData pd = bgy3d_problem_data ();

  /* Petsc does real work: */
  bgy3d_snes_default (&pd, ctx, F, x);

  vec_destroy (&x);       /*  should  not free() */
}
