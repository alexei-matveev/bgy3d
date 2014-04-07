/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

#include <assert.h>

#include "bgy3d.h"
#include "bgy3d-vec.h"          /* vec_destroy() */
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-fft.h"
#include <complex.h>            /* after fftw.h */

double bgy3d_fft_test (int m, int n, int p)
{
  const int N[3] = {m, n, p};
  const int NNN = N[0] * N[1] * N[2];
  Mat A;
  DA da, dc;

  bgy3d_fft_mat_create (N, &A, &da, &dc);

  Vec x = vec_create (da); /* real */
  Vec z = vec_create (da); /* real */

  /* This one is complex, note use of another array descriptor: */
  Vec y = vec_create (dc); /* complex */

  /*
    To  test  if  reference  counting  works, let  us  destroy  array
    descriptors right away. They are  still referenced by and used by
    FFT matrix:
  */
  DMDestroy (&da);
  DMDestroy (&dc);

  VecSetRandom (x, NULL);
  /* VecSet (x, 1.0); */

  /* This corresponds to direct FFT, y = fft(x): */
  MatMult (A, x, y);

  /* This corresponds to inverse FFT, z = ifft(y): */
  MatMultTranspose (A, y, z);

  /*
    The matrix may  be made orthogonal so that  the intuitive relation
    A^T A * x ==  x holds. At the moment A^T * A =  V with V being the
    grid "volume" (number of points):
  */
  VecScale (z, 1.0 / NNN);
  /* VecView (z, PETSC_VIEWER_STDOUT_WORLD); */
  VecAXPY (z, -1.0, x);

  double norm;
  VecNorm (z, NORM_INFINITY, &norm);

  vec_destroy (&x);
  vec_destroy (&y);
  vec_destroy (&z);

  /*
    Also   destroys   array   descriptors   for  real   and   complex
    vectors.  With proper reference  counting the  finalization order
    should not matter:
  */
  MatDestroy (&A);               /* FIXME: petsc 3.2 takes the address! */

  return norm;
}
