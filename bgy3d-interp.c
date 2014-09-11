/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2014 Alexei Matveev

  Trilinear interpolation, eventually cheaper and free of ringing when
  compared to the trigonometric interpolation in bgy3d-fftw3.c
*/

#include "bgy3d.h"              /* comm_allreduce() and such */
#include <assert.h>             /* FIXME: remove asserts() */
#include "bgy3d-vec.h"          /* da_ref() */
#include "bgy3d-interp.h"


/*
  The   result  of   operator  %   may  be   negative  (implementation
  defined). So that  C needs a separate modulo  operation. But it sure
  holds that (a / b) * b + (a % b) = a.  This one will not work as you
  may expect for negative b:
*/
static inline int
mod (int a, int b)
{
  return (a % b + b) % b;
}


/*
  Given  the  real Vec  V  of some  equally-spaced  grid  sample of  a
  real-valued function  v(r), compute the values of  v(r) at arbitrary
  given points  r[][] by trilinear interpolation. Note  that the scale
  of x  here is  such that  the components range  between 0.0  and the
  respective N[]  within the unit cell.   FIXME: Extrapolation outside
  the cell is ...
*/
void
bgy3d_interp (DA da, const Vec V,      /* real, in */
              int np, double r[np][3], /* in */
              double v[np])            /* out */
{
  /* Global dimensions: */
  int N[3];
  da_shape (da, N);

  /* Know the local portion of grid */
  int c[3], n[3];             /* corner and the size of the section */
  DMDAGetCorners (da, &c[0], &c[1], &c[2], &n[0], &n[1], &n[2]);

  real ***V_;
  DMDAVecGetArray (da, V, &V_);

  /* Return v(i, j, k) stored at  V_[k][j][i], or zero if out of local
     range: */
  inline real f (int i, int j, int k)
  {
    /* Periodic extrapolation here: */
    i = mod (i, N[0]);
    j = mod (j, N[1]);
    k = mod (k, N[2]);

    /* FIXME: These are really paranoid now --- redundant if you trust
       mod(): */
    assert (0 <= i && i < N[0]);
    assert (0 <= j && j < N[1]);
    assert (0 <= k && k < N[2]);

    const bool ok =
      c[0] <= i && i < c[0] + n[0] &&
      c[1] <= j && j < c[1] + n[1] &&
      c[2] <= k && k < c[2] + n[2];

    /* Note the sequence, [k][j][i]: */
    return ok? V_[k][j][i] : 0.0;
  }

  for (int p = 0; p < np; p++)
    {
      /*
        Note that even  with the periodic extrapolation in  f(i, j, k)
        the  case of  negative coordinates  will have  a gotcha,  as a
        simple  casting of  a double  to  integer does  not produce  a
        positive remainder, because of rounding toward zero.

        Round  them down instead.  Here (i,  j, k)  is the  left lower
        corner of the volume element  the point p falls into.  And (x,
        y,  z)  are corresponding  fractional  coordinates within  the
        volume element, which shoud be between 0 and 1.

        In an unlikely case when e.g. -1.9 was cast to i = -1 and thus
        x = -0.9  we turn them into i  = -2 and x = 0.1.  That is what
        floor() is for.
      */
      const int i = floor (r[p][0]);
      const int j = floor (r[p][1]);
      const int k = floor (r[p][2]);

      const real x = r[p][0] - i;
      const real y = r[p][1] - j;
      const real z = r[p][2] - k;

      /* Interpolate four corners of a yz-plane section along x: */
      const real c00 =
        f (i + 0, j + 0, k + 0) * (1 - x) +
        f (i + 1, j + 0, k + 0) * x;

      const real c10 =
        f (i + 0, j + 1, k + 0) * (1 - x) +
        f (i + 1, j + 1, k + 0) * x;

      const real c01 =
        f (i + 0, j + 0, k + 1) * (1 - x) +
        f (i + 1, j + 0, k + 1) * x;

      const real c11 =
        f (i + 0, j + 1, k + 1) * (1 - x) +
        f (i + 1, j + 1, k + 1) * x;

      /* Interpolate two ends of a line section along y: */
      const real c0 = c00 * (1 - y) + c10 * y;
      const real c1 = c01 * (1 - y) + c11 * y;

      /* Interpolate a point along z */
      const real c = c0 * (1 - z) + c1 * z;

      v[p] = c;
    }

  /* Used while computing f(i, j, k): */
  DMDAVecRestoreArray (da, V, &V_);

  /*
    Note that  the 8 values  f(i, j, k)  at the corners of  the volume
    element enter the final expression  in linear fashion.  So if each
    worker when not having a  value for some corners uses 0.0 instead,
    than  after adding  terms from  all  workers the  values are  fine
    again:
  */
  comm_allreduce (np, v);
}
