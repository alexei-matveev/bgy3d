/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */

#ifdef WITH_EXTRA_SOLVERS
/* Deprecated.   Used only  in  the older  solvers.   Relies on  fft3d
   implementation in ./fft: */
FFT_DATA *ComputeFFTfromVec (DA da, struct fft_plan_3d *fft_plan, Vec g,
                             FFT_DATA *g_fft);
void ComputeVecfromFFT (DA da, struct fft_plan_3d *fft_plan, Vec g,
                        FFT_DATA *g_fft);
#endif

double bgy3d_fft_test (int m, int n, int k);
