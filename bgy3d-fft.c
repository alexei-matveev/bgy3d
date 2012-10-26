/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.15 2007-04-23 16:55:13 jager Exp $ */
/*==========================================================*/

#include <assert.h>
#include <fftw_mpi.h>
#include "fft_3d.h"             /* FFT_DATA */
#include "petscda.h"            /* DA, Vec */
#include "bgy3d-fftw.h"
#include "bgy3d-fft.h"

#ifdef WITH_EXTRA_SOLVERS
static void unpack (DA da, Vec g, fftw_complex *restrict g_fft)
{
  int index, i0, j0, k0, ni, nj, nk;
  PetscScalar ***g_vec;

  /* Get local portion of the grid */
  DAGetCorners(da, &i0, &j0, &k0, &ni, &nj, &nk);

  DAVecGetArray(da, g, &g_vec);

  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index = 0;
  for (int k = k0; k < k0 + nk; k++)
    for (int j = j0; j < j0 + nj; j++)
      for (int i = i0; i < i0 + ni; i++)
        {
          g_fft[index].re = g_vec[k][j][i];
          g_fft[index].im = 0;  /* Vec g is real */
          index++;
        }
  DAVecRestoreArray(da, g, &g_vec);
}

static void pack (DA da, Vec g, const fftw_complex *restrict g_fft)
{
  int index, i0, j0, k0, ni, nj, nk;
  PetscScalar ***g_vec;

  /* Get local portion of the grid */
  DAGetCorners(da, &i0, &j0, &k0, &ni, &nj, &nk);

  DAVecGetArray(da, g, &g_vec);

  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index = 0;
  for (int k = k0; k < k0 + nk; k++)
    for (int j = j0; j < j0 + nj; j++)
      for (int i = i0; i < i0 + ni; i++)
        {
          /* FIXME: this is subtle, we are packing complex vector into
             a real array. Imaginary part gets ignored: */
          g_vec[k][j][i] = g_fft[index].re;
          index++;
        }
  DAVecRestoreArray(da, g, &g_vec);
}
#endif

#ifdef WITH_EXTRA_SOLVERS
FFT_DATA *ComputeFFTfromVec(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft)
{
  int x[3], n[3];

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  if(g_fft==NULL)
    g_fft = (FFT_DATA*) calloc(n[0] * n[1] * n[2], sizeof(*g_fft));

  /* Real Vec into complex array: */
  unpack (da, g, g_fft);

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
  pack (da, g, g_fft);
}
#endif


double bgy3d_fft_test (int m, int n, int p)
{
  const int N[3] = {m, n, p};
  const int NNN = N[0] * N[1] * N[2];
  Mat A;
  DA da, dc;

  bgy3d_fft_mat_create (N, &A, &da, &dc);

  Vec x, z;                     /* real */
  Vec y;                        /* complex */

  DACreateGlobalVector (da, &x);
  DACreateGlobalVector (da, &z);

  /* This one is complex, note use of another array descriptor: */
  DACreateGlobalVector (dc, &y);

  VecSetRandom (x, NULL);
  /* VecSet (x, 1.0); */

  /* This corresponds to direct FFT, y = fft(x): */
  MatMult (A, x, y);

  /* This corresponds to inverse FFT, z = fft(y): */
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

  VecDestroy (x);
  VecDestroy (y);
  VecDestroy (z);

  /* FIXME:  Also  destroys array  descriptors  for  real and  complex
     vectors: */
  MatDestroy (A);

  DADestroy (da);
  DADestroy (dc);

  return norm;
}
