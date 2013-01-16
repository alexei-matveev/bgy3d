/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-vec.h"

static int vec_local_size (Vec xs)
{
  int n;
  VecGetLocalSize (xs, &n);

  return n;
}

/* ys  = map  (f, xs).  Should also  work with  aliased  arguments for
   in-place transform: */
void bgy3d_vec_map1 (Vec ys, real (*f)(real x), Vec xs)
{
  real *xs_, *ys_;

  VecGetArray (xs, &xs_);
  VecGetArray (ys, &ys_);

  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);

  for (int i = 0; i < n; i++)
    ys_[i] = f (xs_[i]);

  VecRestoreArray (xs, &xs_);
  VecRestoreArray (ys, &ys_);
}

/* zs = map (f, xs, ys).   Should also work with aliased arguments for
   in-place transform: */
void bgy3d_vec_map2 (Vec zs, real (*f)(real x, real y), Vec xs, Vec ys)
{
    real *xs_, *ys_, *zs_;

  VecGetArray (xs, &xs_);
  VecGetArray (ys, &ys_);
  VecGetArray (zs, &zs_);

  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);
  assert (vec_local_size (zs) == n);

  for (int i = 0; i < n; i++)
    zs_[i] = f (xs_[i], ys_[i]);

  VecRestoreArray (xs, &xs_);
  VecRestoreArray (ys, &ys_);
  VecRestoreArray (zs, &zs_);
}


/* This one is supposed to save enough meta-info (such as distribution
   pattern, dimensions) to recover the vector from scratch: */
void bgy3d_vec_save (const char file[], const Vec vec)
{
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_WRITE, &viewer);
  VecView (vec, viewer);
  PetscViewerDestroy (viewer);
}


/* This one  returns a newly  created vector: */
Vec bgy3d_vec_load (const char file[])
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
void bgy3d_vec_read (const char file[], Vec vec)
{
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
  VecLoadIntoVector (viewer, vec);
  PetscViewerDestroy (viewer);
}


/* This one will not save much of a meta-info: */
void bgy3d_vec_save_ascii (const char file[], const Vec vec)
{
  PetscViewer viewer;

  PetscViewerASCIIOpen (PETSC_COMM_WORLD, file, &viewer);
  /* PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_MATLAB); */
  /* PetscViewerSetFormat (viewer, PETSC_VIEWER_DEFAULT); */
  /* PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_VTK); */
  PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_COMMON);
  VecView (vec, viewer);
  PetscViewerDestroy (viewer);
}
