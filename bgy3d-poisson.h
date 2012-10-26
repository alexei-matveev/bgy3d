/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */

void bgy3d_poisson (const State *BHD, Vec uc, Vec rho, real q);

#ifdef L_BOUNDARY
void InitializeLaplaceMatrix (const State *BHD, real zpad);
void InitializeKSPSolver (State *BHD);
real ImposeLaplaceBoundary (const State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif
