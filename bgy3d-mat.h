/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013 Alexei Matveev
*/

#include "bgy3d.h"

/* Functions that take a Mat and a Vec and update another Vec: */
typedef PetscErrorCode (*Operation)(Mat A, Vec x, Vec y);

/* Functions that take a Mat and destroy its internals: */
typedef PetscErrorCode (*Destructor)(Mat A);


/* Returns an inverse matrix. */
Mat mat_inverse (Mat A);


/* Petsc  3.2 changed the  interface of  XXXDestroy() methods  so that
   they take the pointer to a Petsc object and nullify it: */
static inline
void mat_destroy (Mat *A)
{
  /* Since Petsc 3.2, MatDestroy() also zeroed out the buffer and
   * takes the address instead of vector as input argument */
  MatDestroy (A);
  /* FIXME: only needed before Petsc 3.2 */
  *A = NULL;
}


/*
  Convenience functions  to deal  with matrix free  implementations of
  linear operators.

  Untyped  convinience wrapper.   Not all  matrices have  a meaningful
  context:
*/
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
