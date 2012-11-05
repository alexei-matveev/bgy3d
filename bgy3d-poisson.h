/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */

void Zeropad_Function (const State *BHD, Vec g, real shift);

void bgy3d_poisson (const State *BHD, Vec uc, Vec rho, real q);

#ifdef L_BOUNDARY
void bgy3d_laplace_create (const DA da, const ProblemData *PD, Mat *M, KSP *ksp);
void ImposeLaplaceBoundary (const State *BHD, Vec v, Vec b, Vec x);
#endif
