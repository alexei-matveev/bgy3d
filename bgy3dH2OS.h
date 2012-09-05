/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
void ReadPairDistribution (const State *BHD, const char *filename, Vec g2);
Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini);
Vec BGY3dM_solve_H2O_2site(const ProblemData *PD, Vec g_ini);
void RecomputeInitialFFTs (State *BHD, real damp, real damp_LJ);
void RecomputeInitialSoluteData(State *BHD, real damp, real damp_LJ, real zpad);
real ComputeCharge (State *BHD, Vec g1, Vec g2);
void Compute_H2O_interS (const State *BHD,
                         fftw_complex *(fg2_fft[3]), Vec g, real rho, Vec dg_help);
#ifdef L_BOUNDARY
void InitializeLaplaceMatrix (State *BHD, real zpad);
void InitializeKSPSolver (State *BHD);
real  ImposeLaplaceBoundary (const State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif

