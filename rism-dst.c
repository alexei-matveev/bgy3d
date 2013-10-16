/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include <fftw3.h>
#include <assert.h>
#include "rism-dst.h"

void rism_dst (size_t n, double out[n], const double in[n])
{
  /* FIXME: does it write to in[]? */
  fftw_plan plan = fftw_plan_r2r_1d (n, (double*) in, out, FFTW_RODFT11, FFTW_ESTIMATE);

  fftw_execute (plan);

  fftw_destroy_plan (plan);
}

/* Transform m continous arrays each  of length n. In Fortran terms do
   FFT for each column of the n x m matrix buf(:, :). */
void rism_dst_columns (int m, int n, double buf[m][n])
{
  const fftw_r2r_kind kind = FFTW_RODFT11;

  fftw_plan plan =
    fftw_plan_many_r2r (1, &n, m, /* rank, dimensions[rank], howmany */
                        (double*) buf, NULL, 1, n, /* inp, ?, istride, idist */
                        (double*) buf, NULL, 1, n, /* out, ?, ostride, odist */
                        &kind, FFTW_ESTIMATE); /* kind, flags */
  assert (plan != NULL);

  fftw_execute (plan);

  fftw_destroy_plan (plan);
}


/* Transform m stride-m  arrays each of length n.  In Fortran terms do
   FFT for each row of the m x n matrix buf(:, :). */
void rism_dst_rows (int n, int m, double buf[n][m])
{
  const fftw_r2r_kind kind = FFTW_RODFT11;

  fftw_plan plan =
    fftw_plan_many_r2r (1, &n, m, /* rank, dimensions[rank], howmany */
                        (double*) buf, NULL, m, 1, /* inp, ?, istride, idist */
                        (double*) buf, NULL, m, 1, /* out, ?, ostride, odist */
                        &kind, FFTW_ESTIMATE); /* kind, flags */
  assert (plan != NULL);

  fftw_execute (plan);

  fftw_destroy_plan (plan);
}
