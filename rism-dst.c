/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include <fftw3.h>
#include "rism-dst.h"

void rism_dst (size_t n, double out[n], const double in[n])
{
  /* FIXME: does it write to in[]? */
  fftw_plan plan = fftw_plan_r2r_1d (n, (double*) in, out, FFTW_RODFT11, FFTW_ESTIMATE);

  fftw_execute (plan);

  fftw_destroy_plan (plan);
}


void rism_dst_many (int m, int n, double buf[m][n])
{
  const fftw_r2r_kind kind = FFTW_RODFT11;

  fftw_plan plan = fftw_plan_many_r2r (1, &n, m,
                                       (double*) buf, NULL, 1, n,
                                       (double*) buf, NULL, 1, n,
                                       &kind, FFTW_ESTIMATE);

  fftw_execute (plan);

  fftw_destroy_plan (plan);
}
