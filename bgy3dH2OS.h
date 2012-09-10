/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */


/*
  This function is the main entry  point for the BGY3dM equation for a
  2-site solvent and an arbitrary solute.  The two vectors in

  Vec g[2], intent(out)

  are  initialzed as  global distributed  arrays and  filled  with the
  solvent site  distributions. It is the responsibility  of the caller
  to destroy them when no more needed.
 */
void bgy3d_solve_with_solute (const ProblemData *PD,
                              int n, const Site solute[n],
                              Vec g[2]);

void ReadPairDistribution (const State *BHD, const char *filename, Vec g2);
Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini);
Vec BGY3dM_solve_H2O_2site(const ProblemData *PD, Vec g_ini);
void RecomputeInitialFFTs (State *BHD, real damp, real damp_LJ);
real ComputeCharge (State *BHD, Vec g1, Vec g2);
void Compute_H2O_interS (const State *BHD,
                         fftw_complex *(fg2_fft[3]), Vec g, real rho, Vec dg_help);
#ifdef L_BOUNDARY
void InitializeLaplaceMatrix (State *BHD, real zpad);
void InitializeKSPSolver (State *BHD);
real  ImposeLaplaceBoundary (const State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter);
#endif

