/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.15 2007-04-23 16:55:13 jager Exp $ */
/*==========================================================*/

#include <assert.h>

#ifdef WITH_EXTRA_SOLVERS
#include <fftw_mpi.h>           /* to get FFT_DATA set */
#include "fft_3d.h"             /* FFT_DATA */
#endif

#include "bgy3d.h"
#include "bgy3d-vec.h"          /* bgy3d_vec_destroy() */
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-fft.h"
#include <complex.h>            /* after fftw.h */

#ifdef WITH_EXTRA_SOLVERS
static void unpack (DA da, Vec g, complex *restrict g_fft)
{
  /* Get local portion of the grid */
  int i0, j0, k0, ni, nj, nk;
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  PetscScalar ***g_;
  DAVecGetArray (da, g, &g_);

  /* loop over local portion of grid */
  int ijk = 0;
  for (int k = k0; k < k0 + nk; k++)
    for (int j = j0; j < j0 + nj; j++)
      for (int i = i0; i < i0 + ni; i++)
        g_fft[ijk++] = g_[k][j][i]; /* Vec g is real */

  DAVecRestoreArray (da, g, &g_);
}

/* NOTE: this  is subtle,  we are packing  complex vector into  a real
   array. Imaginary part gets ignored: */
static void pack (DA da, Vec g, const complex *restrict g_fft)
{
  /* Get local portion of the grid */
  int i0, j0, k0, ni, nj, nk;
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  PetscScalar ***g_;
  DAVecGetArray (da, g, &g_);

  /* loop over local portion of grid */
  int ijk = 0;
  for (int k = k0; k < k0 + nk; k++)
    for (int j = j0; j < j0 + nj; j++)
      for (int i = i0; i < i0 + ni; i++)
        g_[k][j][i] = g_fft[ijk++]; /* drops imaginary part */

  DAVecRestoreArray (da, g, &g_);
}

FFT_DATA *ComputeFFTfromVec(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft)
{
  int x[3], n[3];

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  if(g_fft==NULL)
    g_fft = (FFT_DATA*) calloc(n[0] * n[1] * n[2], sizeof(*g_fft));

  /* Real Vec into complex array: */
  unpack (da, g, (complex*) g_fft);

  /* forward fft */
  fft_3d(g_fft, g_fft, 1, fft_plan);

  return g_fft;
}


void ComputeVecfromFFT(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft)
{
  assert (g_fft != NULL);

  /* backward fft */
  fft_3d(g_fft, g_fft, -1, fft_plan);

  /* Pack a  complex vector with (hopefully)  vanishing imaginary part
     into a real Vec: */
  pack (da, g, (complex*) g_fft);
}
#endif


double bgy3d_fft_test (int m, int n, int p)
{
  const int N[3] = {m, n, p};
  const int NNN = N[0] * N[1] * N[2];
  Mat A;
  DA da, dc;

  bgy3d_fft_mat_create (N, &A, &da, &dc);

  Vec x = bgy3d_vec_create (da); /* real */
  Vec z = bgy3d_vec_create (da); /* real */

  /* This one is complex, note use of another array descriptor: */
  Vec y = bgy3d_vec_create (dc); /* complex */

  /*
    To  test  if  reference  counting  works, let  us  destroy  array
    descriptors right away. They are  still referenced by and used by
    FFT matrix:
  */
  DADestroy (da);
  DADestroy (dc);

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

  bgy3d_vec_destroy (&x);
  bgy3d_vec_destroy (&y);
  bgy3d_vec_destroy (&z);

  /*
    Also   destroys   array   descriptors   for  real   and   complex
    vectors.  With proper reference  counting the  finalization order
    should not matter:
  */
  MatDestroy (&A);               /* FIXME: petsc 3.2 takes the address! */

  return norm;
}
