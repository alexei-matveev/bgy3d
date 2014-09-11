/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013 Alexei Matveev
*/

#include "bgy3d.h"              /* State */
#include "bgy3d-vec.h"          /* vec_create() */
#include "bgy3d-fftw.h"         /* bgy3d_fft_interp() */
#include "bgy3d-interp.h"       /* bgy3d_interp() */
#include "lebed/lebed.h"        /* genpts() */


/* Flip this to use expensive trigonometric interpolation: */
static bool trilinear = true;

void
rism_rdf (State *dom, Vec g,
          const real a[3],        /* center */
          int n, const real r[n], /* radial mesh */
          int m,                  /* quadrature order */
          real rdf[n])            /* out */
{
  const ProblemData *PD = dom->PD;

  /* Prepare Fourier coefficients.  Only needed for trigonometric
     interpolation: */
  local Vec g_fft = NULL;
  if (!trilinear)
    {
      g_fft = vec_create (dom->dc);

      MatMult (dom->fft_mat, g, g_fft);
      /* Do not VecScale (g_fft, volume_element (PD)), the
         interpolation code assumes that ... */
    }

  /* Coordinates  and weights of  spherical quadrature.  Only mm  <= m
     points are valid: */
  real x[m][3], w[m];

  const int mm = genpts (m, x, w);

  /* (Real) grid coordinates and interpolated values: */
  real y[mm][3], gr[mm];

  /* Over radial layers: */
  for (int j = 0; j < n; j++)
    {
      /* Translate quadrature  coordinates into real  grid coordinates
         where integer values correspond to the grid nodes: */
      for (int i = 0; i < mm; i++)
        FOR_DIM
          y[i][dim] = (a[dim] + r[j] * x[i][dim] + PD->L[dim] / 2) / PD->h[dim];

      /* NOTE:  trigonometric  interpolation   is  O(N^3)  per  point,
         trilinear is O(1): */
      if (trilinear)
        bgy3d_interp (dom->da, g, mm, y, gr);
      else
        bgy3d_fft_interp (dom->fft_mat, g_fft, mm, y, gr);

      /* Scale by integration weights: */
      for (int i = 0; i < mm; i++)
        gr[i] *= w[i];

      /* Integrate over sphere: */
      rdf[j] = sum (mm, gr);
    }

  if (!trilinear)
    vec_destroy (&g_fft);
}
