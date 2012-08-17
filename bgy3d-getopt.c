#include <assert.h>
#include "petscsys.h"
#include "petscda.h" // Vec
#include "bgy3d-getopt.h"

int bgy3d_getopt_test (const char key[])
{
  PetscTruth test;

  /* Petsc insists  on keys  having a leading  dash. The  prefix (here
     PETSC_NULL) cannot have a leading dash. Go figure ... */
  PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULL, key, &test);
  assert (!ierr);

  return (int) test;
}

int bgy3d_getopt_int (const char key[], int *val)
{
  PetscTruth test;

  PetscErrorCode ierr = PetscOptionsGetInt(PETSC_NULL, key, val, &test);
  assert (!ierr);

  return (int) test;
}

int bgy3d_getopt_real (const char key[], double *val)
{
  PetscTruth test;

  PetscErrorCode ierr = PetscOptionsGetReal(PETSC_NULL, key, val, &test);
  assert (!ierr);

  return (int) test;
}


void bgy3d_load_vec (const char file[], Vec *vec) {
    PetscViewer viewer;

    PetscViewerBinaryOpen(PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
    VecLoad(viewer, VECMPI, vec);
    PetscViewerDestroy(viewer);
}
