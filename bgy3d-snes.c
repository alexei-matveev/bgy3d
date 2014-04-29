/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

#include "petscsnes.h"          /* SNESSolve, etc. Before bgy3d.h */
#include "bgy3d.h"
#include "bgy3d-vec.h"          /* vec_duplicate() */
#include "bgy3d-getopt.h"       /* bgy3d_getopt_string() */
#include "bgy3d-mat.h"          /* mat_shell_create() */
#include "bgy3d-snes.h"         /* VecFunc1, ArrFunc1 */


void bgy3d_snes_default (const ProblemData *PD, void *ctx,
                         VecFunc1 F, VecFunc2 dF, Vec x)
{
  char solver[20] = "newton";
  bgy3d_getopt_string ("--snes-solver", sizeof solver, solver);

  if (strcmp (solver, "newton") == 0)
    bgy3d_snes_newton (PD, ctx, F, dF, x);
  else if (strcmp (solver, "picard") == 0)
    bgy3d_snes_picard (PD, ctx, F, dF, x);
  else if (strcmp (solver, "jager") == 0)
    bgy3d_snes_jager (PD, ctx, F, dF, x);
  else if (strcmp (solver, "trial") == 0)
    bgy3d_snes_trial (PD, ctx, F, dF, x);
  else
    {
      PetscPrintf (PETSC_COMM_WORLD, "No such SNES solver: %s\n", solver);
      exit (1);
    }
}


/* Matrix-free matrix  based on linear  VecFunc1. Note it is  the user
   responsibility to guarantee that VecFunc1 is linear. */
typedef struct Mtx
{
  void *data;
  VecFunc1 f;
} Mtx;


/* [Operation] Apply  J. An Operation that  takes a Mat and  a Vec and
   fills another Vec: */
static PetscErrorCode
mapply (Mat J, Vec x, Vec y)
{
  const Mtx *ctx = mat_shell_context (J);

  /* y = J * x */
  ctx->f (ctx->data, x, y);
  return 0;
}


/* [Destructor]   A  Destructor   that  takes   a  Mat   and  destroys
   internals: */
static PetscErrorCode
mdestroy (Mat J)
{
  free (mat_shell_context (J)); /* malloc() in mcreate() */
  return 0;
}


/* [Constructor] */
static Mat
mcreate (int n, void *data, VecFunc1 f)
{
  Mtx *ctx = malloc (sizeof *ctx); /* free() in mdestroy() */

  *ctx = (struct Mtx) {data, f};

  /* A square matrix with one  Operation and Destructor. Here n is the
     size of the local vector section the matrix operates upon. */
  return mat_shell_create (n, ctx, mapply, mdestroy);
}


/* Matrix-free  Jacobian,   J(r).   The  content  of  the   Vec  r  is
   occasionally updated by PETSC invoking jupdate(). */
typedef struct Ctx
{
  void *data;
  Vec r;
  VecFunc2 j;
} Ctx;


/* [Operation] Apply J(r): */
static PetscErrorCode
japply (Mat J, Vec x, Vec y)
{
  const Ctx *ctx = mat_shell_context (J);

  /* y = J(r) * x */
  ctx->j (ctx->data, ctx->r, x, y);
  return 0;
}


/* [Destructor] Destroy internals of J(r): */
static PetscErrorCode
jdestroy (Mat J)
{
  Ctx *ctx = mat_shell_context (J);
  vec_destroy (&ctx->r);        /* vec_ref() in jcreate() */
  free (ctx);                   /* malloc() in jcreate() */
  return 0;
}


/* [Constructor] */
static Mat
jcreate (void *data, Vec r, VecFunc2 j)
{
  Ctx *ctx = malloc (sizeof *ctx); /* free() in jdestroy() */

  /* Inrement  referenc count,  to  allow the  caller  to destroy  the
     Vec. The corresponding vec_destroy() in jdestroy(): */
  *ctx = (struct Ctx) {data, vec_ref (r), j};

  /* A square matrix with one Operation and Destructor: */
  return mat_shell_create (vec_local_size (r), ctx, japply, jdestroy);
}


/*
  This  one is  called  by the  SNES  solver to  update (re-seat)  the
  Jacobian   J(r)   ->   J(r').    The   interface   is   rigid,   see
  MatMFFDComputeJacobian()
*/
static PetscErrorCode
jupdate (SNES snes, Vec r, Mat *J, Mat *B, MatStructure *s, void *data)
{
  (void) snes;
  (void) B;
  (void) s;
  (void) data;

  Ctx *ctx = mat_shell_context (*J);
  VecCopy (r, ctx->r);

  return 0;
}



/*
  For solving HNC equation with Newton. Except of x everything else in
  the  closure context is  considered read  only input  or intemediate
  terms depending  on x.  That  is when looking for  total correlation
  function h, the direct correlation function should be fixed (or be a
  function of h).
*/
void
bgy3d_snes_newton (const ProblemData *PD, void *ctx,
                   VecFunc1 F, VecFunc2 dF, Vec x)
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

  /* Line search: SNESNEWTONLS (aka SNESLS) trust region: SNESTR */
  SNESSetType (snes, SNESNEWTONLS);

  if (dF)
    {
      /* The point  to compute (rather  apply) Jacobian will  be saved
         here: */
      local Vec r0 = vec_duplicate (r);

      /* Self-made  matrix-free Jacobian  based  on the  user-supplied
         VecFunc2:  */
      local Mat J = jcreate (ctx, r0, dF);

      /* Stuck it into SNES object: */
      SNESSetJacobian (snes, J, J, jupdate, ctx);

      /* Both  Mat  J  and  SNES  snes  should  have  incremented  the
         refcounts: */
      mat_destroy (&J);
      vec_destroy (&r0);
    }
  else
    {
      /*
        This  has the  same  effect as  specifying  "-snes_mf" in  the
        command  line.  We  do  it  explicitly here  in  order not  to
        require that switch from the user. MF stays for Matrix-Free.

        The  idea is  that  matrix-vector products  with Jacobian  are
        computed   by  numerical   differentiation  of   the  original
        form-function.  The form-function has to already be associated
        with the SNES object, see above.
      */
      local Mat J;              /* I guess we need to destroy it? */
      MatCreateSNESMF (snes, &J);

      /*
        Petsc  provides  a convenience  function  for Jacobian  update
        which cooperates  with the  matrix-free Jacobian J  created by
        MatCreateSNESMF() as above:

        PetscErrorCode
        MatMFFDComputeJacobian (SNES, Vec, Mat*, Mat*, MatStructure*, void*)

        The last  argument is a user context  for Jacobian evaluation.
        This function is  only used to *update* the  Jacobian.  In the
        matrix-free  case an  update amounts  to noting  (copying) the
        location Vec r  at which the Jacobian J(r)  is to be evaluated
        next  time it  is applied  as in  J(r) *  dr.  The application
        itself is performed by Mat J matrix shell.
      */
      SNESSetJacobian (snes, J, J, MatMFFDComputeJacobian, NULL);
      mat_destroy (&J);     /* I hope SNES saved a ref to that Mat? */
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

    /* Set rtol, atol, dtol, maxits */
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
        here (compare to the  motivation and settings in bgy3d-mat.c).
        With the  defaults above the  code below requires the  norm of
        the  residual  of  the  *linear* Jacobian  sub-problem  to  be
        reduced by at least factor 1000  in order to be used as search
        direction:
      */
      KSPSetTolerances (ksp, 100 * rtol, abstol, 10 * dtol, maxits / 10);

      /*
        A call to  KSPSetFromOptions() will eventually overwrite these
        settigns   by  the   runtime  options   (including   those  in
        ~/.petscrc).  Presumably,  an equivalent of this  call is made
        from the body of SNESSetFromOptions(), see below.
      */
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

  /* FIXME: what  should we do  if there is no  convergence?  Negative
     value indicates diverged, positive value converged: */
  SNESConvergedReason reason;
  SNESGetConvergedReason (snes, &reason);
  SNESDestroy (&snes);

  if (reason <= 0)
    misc_error (__func__, SNESConvergedReasons[reason]); /* longjmp! */
  assert (reason > 0);
}

void bgy3d_snes_picard (const ProblemData *PD, void *ctx,
                        VecFunc1 F, VecFunc2 dF, Vec x)
{
  if (dF)
    fprintf (stderr, "bgy3d_snes_picard: Warning: not using Jacobian!");

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


void bgy3d_snes_trial (const ProblemData *PD, void *ctx,
                       VecFunc1 F, VecFunc2 dF, Vec x)
{
  /* First do a few "slow" iterations: */
  {
    ProblemData pd = *PD;       /* modify a copy, not the original */
    pd.max_iter = 100;
    pd.lambda = 0.02;
    bgy3d_snes_picard (&pd, ctx, F, dF, x);
  }
  /* Then continue with Newton: */
  bgy3d_snes_newton (PD, ctx, F, dF, x);
}


void bgy3d_snes_jager (const ProblemData *PD, void *ctx,
                       VecFunc1 F, VecFunc2 dF, Vec x)
{
  if (dF)
    fprintf (stderr, "bgy3d_snes_jager: Warning: not using Jacobian!");

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
  A solver  for an untyped  (array) form-function ArrFunc1.   Finds an
  x_[n] such that  dx_[n] as returned by ArrFunc1 f  (ctx, n, x_, dx_)
  is zero.  May  be eventually used from Fortran,  so if ever changing
  the interface  update the  corresponding interface block  in Fortran
  sources.

  FIXME:  in parallel  runs with  P workers  this function  appears to
  operate  with a  distributed Vec  of  total length  n *  P built  of
  redundant sections of length n on each worker.
*/
void rism_snes (void *ctx, ArrFunc1 f, ArrFunc2 df, int n, real x_[n])
{
  local Vec x = vec_from_array (n, x_);

  /* Implements VecFunc1 interface: */
  void F (void *ctx, Vec y, Vec dy)
  {
    local real *y_ = vec_get_array (y);
    local real *dy_ = vec_get_array (dy);

    /* Implements ArrFunc1 interface: */
    f (ctx, n, y_, dy_);

    vec_restore_array (y, &y_);
    vec_restore_array (dy, &dy_);
  }

  /* Implements VecFunc2 interface: */
  void dF (void *ctx, Vec a, Vec b, Vec c)
  {
    local real *a_ = vec_get_array (a);
    local real *b_ = vec_get_array (b);
    local real *c_ = vec_get_array (c);

    /* Implements ArrFunc2 interface: */
    df (ctx, n, a_, b_, c_);

    vec_restore_array (a, &a_);
    vec_restore_array (b, &b_);
    vec_restore_array (c, &c_);
  }

  /* The result depends on the command line and affects the solver: */
  ProblemData pd = bgy3d_problem_data ();

  /* Petsc does real work: */
  if (df)
    bgy3d_snes_default (&pd, ctx, F, dF, x);
  else
    bgy3d_snes_default (&pd, ctx, F, NULL, x); /* dF is never NULL */

  vec_destroy (&x);       /*  should  not free() */
}


/* For solving linear equation F x = b iteratively. */
static void
krylov (void *ctx, VecFunc1 F, Vec b, Vec x)
{
  local Mat A = mcreate (vec_local_size (b), ctx, F);

  /* B := A^-1 */
  local Mat B = mat_inverse (A);
  mat_destroy (&A);

  /* x := B * b == A^-1 * b */
  MatMult (B, b, x);

  mat_destroy (&B);
}


/*
  Assumes f(x) is linear and solves for f(x) = b.

  FIXME: adapt for multiple RHSs.

  FIXME:  in parallel  runs with  P workers  this function  appears to
  operate  with a  distributed Vec  of  total length  n *  P built  of
  redundant sections of length n on each worker.
*/
void
rism_krylov (void *ctx, ArrFunc1 f, int n, real b_[n], real x_[n])
{
  local Vec b = vec_from_array (n, b_);
  local Vec x = vec_from_array (n, x_);

  /* Implements VecFunc1 interface: */
  void F (void *ctx, Vec y, Vec dy)
  {
    local real *y_ = vec_get_array (y);
    local real *dy_ = vec_get_array (dy);

    /* Implements ArrFunc1 interface: */
    f (ctx, n, y_, dy_);

    vec_restore_array (y, &y_);
    vec_restore_array (dy, &dy_);
  }

  /* Petsc does the real work: */
  krylov (ctx, F, b, x);

  /* Should not free() */
  vec_destroy (&x);
  vec_destroy (&b);
}
