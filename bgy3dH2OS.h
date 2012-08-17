void ReadPairDistribution (const State *BHD, const char *filename, Vec g2);
Vec BGY3dM_solve_H2O_3site(ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3dM_solve_H2O_2site(ProblemData *PD, Vec g_ini, int vdim);
void RecomputeInitialFFTs (State *BHD, real damp, real damp_LJ);
void RecomputeInitialSoluteData(State *BHD, real damp, real damp_LJ, real zpad);
real ComputeCharge (State *BHD, Vec g1, Vec g2);
void Compute_H2O_interS (State *BHD,
                       fftw_complex *(fg2_fft[3]), Vec g, fftw_complex *coul_fft,
                       fftw_complex *(fs_fft[3]), real con, real rho, Vec dg_help);
#ifdef L_BOUNDARY
void InitializeLaplaceMatrix (State *BHD, real zpad);
void InitializeKSPSolver (State *BHD);
real  ImposeLaplaceBoundary (const State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif

