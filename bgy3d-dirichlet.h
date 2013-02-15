/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#ifdef L_BOUNDARY
Mat bgy3d_dirichlet_create (const DA da, const ProblemData *PD);
void bgy3d_impose_laplace_boundary (const State *BHD, Vec v, Vec b, Vec x);
#endif

Mat bgy3d_lap_mat_create (const DA da, const real h[3]);

