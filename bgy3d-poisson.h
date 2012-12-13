/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void bgy3d_poisson (const State *BHD, Vec uc, Vec rho, real q);

#ifdef L_BOUNDARY
void bgy3d_laplace_create (const DA da, const ProblemData *PD, Mat *M, KSP *ksp);
void bgy3d_impose_laplace_boundary (const State *BHD, Vec v, Vec b, Vec x);
#endif

void bgy3d_lap_mat_create (const DA da, const real h[3], Mat *A);

