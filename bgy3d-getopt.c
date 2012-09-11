/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
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

int bgy3d_getopt_string (const char key[], char *val, size_t len)
{
  PetscTruth test;

  PetscErrorCode ierr = PetscOptionsGetString(PETSC_NULL, key, val, len, &test);
  assert (!ierr);

  return (int) test;
}

/* This one is supposed to save enough meta-info (such as distribution
   pattern, dimensions) to recover the vector from scratch: */
void bgy3d_save_vec (const char file[], const Vec vec)
{
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_WRITE, &viewer);
  VecView (vec, viewer);
  PetscViewerDestroy (viewer);
}

/* This one  returns a newly  created vector: */
Vec bgy3d_load_vec (const char file[])
{
  Vec vec;                      /* new one */
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
  VecLoad (viewer, VECMPI, &vec); /* creates it */
  PetscViewerDestroy (viewer);

  return vec;
}

/* This one takes an allocated vector  and fills it with the data read
   from from disk. Pass a valid vector here: */
void bgy3d_read_vec (const char file[], Vec vec)
{
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
  VecLoadIntoVector (viewer, vec);
  PetscViewerDestroy (viewer);
}

/* This one will not save much of a meta-info: */
void bgy3d_save_vec_ascii (const char file[], const Vec vec)
{
  PetscViewer viewer;

  PetscViewerASCIIOpen (PETSC_COMM_WORLD, file, &viewer);
  PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_MATLAB);
  /* PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_VTK); */
  VecView (vec, viewer);
  PetscViewerDestroy (viewer);
}
