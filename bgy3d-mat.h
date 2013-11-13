/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
#include "bgy3d.h"

/*
  Convenience   types   and    functions   to   create   matrix   free
  implementations of linear operators:
*/

/* Functions that take a Mat and a Vec and update another Vec: */
typedef PetscErrorCode (*Operation)(Mat A, Vec x, Vec y);

/* Functions that take a Mat and destroy its internals: */
typedef PetscErrorCode (*Destructor)(Mat A);


/* Untyped convinience  wrapper.  Not  all matrices have  a meaningful
   context: */
static inline void*
mat_shell_context (Mat A)
{
  void *ctx;
  MatShellGetContext (A, &ctx);
  return ctx;
}


static inline Mat
mat_shell_create (int n, void *ctx, Operation mult, Destructor destroy)
{
  /* Create  matrix shell  with  proper dimensions  and associate  the
     context with it: */
  Mat A;
  MatCreateShell (PETSC_COMM_WORLD,
                  n, n, PETSC_DETERMINE, PETSC_DETERMINE,
                  ctx, &A);

  /* Set matrix operations: */
  MatShellSetOperation (A, MATOP_MULT, (void (*)(void)) mult);
  MatShellSetOperation (A, MATOP_DESTROY, (void (*)(void)) destroy);

  return A;
}
