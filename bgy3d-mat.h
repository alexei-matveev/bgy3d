/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
#include "bgy3d.h"

/*
  Convenience   types   and    functions   to   create   matrix   free
  implementations of linear operators:
 */

/* Functions that take a Mat and a Vec and update another Vec: */
typedef PetscErrorCode (*Operation)(Mat A, Vec x, Vec y);

/* Functions that take a Mat and destroy them: */
typedef PetscErrorCode (*Destructor)(Mat A);


/* Untyped  convinience wrapper,  not all  matrices have  a meaningful
   context: */
static void* context (Mat A)
{
  void *ctx;
  MatShellGetContext (A, &ctx);
  return ctx;
}


/* Creates a matix  shell, but does not associate  any operations with
   it: */
static Mat mat_create_shell (int n, int N, void *ctx)
{
  Mat A;
  MatCreateShell (PETSC_COMM_WORLD, n, n, N, N, ctx, &A);
  return A;
}


static Mat mat_create (int n, int N, void *ctx,
                       Operation mat_mult,
                       Destructor mat_destroy)
{
  /* Create  matrix shell  with  proper dimensions  and associate  the
     context with it: */
  Mat A = mat_create_shell (n, N, ctx);

  /* Set matrix operations: */
  MatShellSetOperation (A, MATOP_MULT,
                        (void (*)(void)) mat_mult);

  MatShellSetOperation (A, MATOP_DESTROY,
                        (void (*)(void)) mat_destroy);

  return A;
}
