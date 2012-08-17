void ReadPairDistribution(BGY3dH2OData BHD, char *filename, Vec g2);
Vec BGY3dM_solve_H2O_3site(ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3dM_solve_H2O_2site(ProblemData *PD, Vec g_ini, int vdim);
void RecomputeInitialFFTs(BGY3dH2OData BHD, real damp, real damp_LJ);
void RecomputeInitialSoluteData(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad);
real ComputeCharge(BGY3dH2OData BHD, Vec g1, Vec g2);
void Compute_H2O_interS(BGY3dH2OData BHD,
                       fftw_complex *(fg2_fft[3]), Vec g, fftw_complex *coul_fft,
                       fftw_complex *(fs_fft[3]), real con, real rho, Vec dg_help);
#ifdef L_BOUNDARY
void InitializeLaplaceMatrix(BGY3dH2OData BHD, real zpad);
void InitializeKSPSolver(BGY3dH2OData BHD);
real  ImposeLaplaceBoundary(BGY3dH2OData BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif

