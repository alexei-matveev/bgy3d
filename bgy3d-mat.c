/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-mat.h"


/*
  To create a matrix one needs  the local and total size of the matrix
  (section).  This  is how  to get it  from another matrix.   See also
  asizes() in bgy3d-dirichlet.c.
*/
static void
msizes (const Mat A, int *n3, int *N3)
{
  int m, n;
  MatGetLocalSize (A, &m, &n);
  assert (m == n);

  int M, N;
  MatGetSize (A, &M, &N);
  assert (M == N);

  *n3 = n;
  *N3 = N;
}


/*
 * Implementation of the matrix inverse by an iterative solver.
 */

/* KSP stays  for Krylov  Sub-Space, used to  solve systems  of linear
   equations: */
static KSP
ksp_create (Mat M)
{
  PC pc;
  KSP ksp;

  /* Create ksp environment */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  KSPGetPC (ksp, &pc);

  /* FIXME: literal tolerances here: */
  KSPSetTolerances (ksp, 1.0e-4, 1.0e-4, 1.0e+5, 1000);

  /* Set Matrix */
  //KSPSetOperators (*ksp, M, M, SAME_NONZERO_PATTERN);
  KSPSetOperators (ksp, M, M, SAME_PRECONDITIONER);

  /* Set preconditioner */
  PCSetType (pc, PCBJACOBI);

  /* This is the place which is  supposed to tell KSP solver to re-use
     the supplied vector as initial guess: */
  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  /* Runtime  options (including  those in  ~/.petscrc)  will override
     default parameters: */
  KSPSetFromOptions (ksp);

  return ksp;
}


/*
  [Operation] Does y =  A x where A = B^-1 and  was constructed as A =
  mat_inverse (B):
*/
static PetscErrorCode
mat_mult_inv (Mat A, Vec x, Vec y)
{
  /* KSP solver was created on matrix construction time: */
  KSP ksp = mat_shell_context (A);

  KSPSolve (ksp, x, y);

  /* If you think you need to know how many iteration were required to
     solve the equation do this: */
  if (verbosity > 0)
    {
      int iter;
      KSPGetIterationNumber (ksp, &iter);
      PetscPrintf (PETSC_COMM_WORLD, "ksp(%2d) ", iter);
    }

  return 0;
}


/*
  [Destructor]   Destroys   internals   of   a   matrix   created   by
  mat_inverse().
*/
static PetscErrorCode
mat_destroy_inv (Mat A)
{
  if (verbosity > 0)
    printf ("mat_destroy_inv(%p)\n", A);

  KSP ksp = mat_shell_context (A);

  KSPDestroy (&ksp);

  return 0;
}


/* [Constructor] Retuns an inverse matrix. */
Mat
mat_inverse (Mat A)
{
  /* An  inverse operator  needs to  solve linear  equations  with the
     original matrix.  Presumable KSP object will hold a reference: */
  KSP ksp = ksp_create (A);

  /* Get the shape of the future matrix: */
  int n, N;
  msizes (A, &n, &N);

  return mat_shell_create (n, ksp, mat_mult_inv, mat_destroy_inv);
}
