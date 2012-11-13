/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */


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
                              void (*density)(int k, const real x[k][3], real rho[k]),
                              Vec g[2],
                              Context **v);

void ReadPairDistribution (const State *BHD, const char *filename, Vec g2);
Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini);
Vec BGY3dM_solve_H2O_2site(const ProblemData *PD, Vec g_ini);
void RecomputeInitialFFTs (State *BHD, real damp, real damp_LJ);
void Compute_H2O_interS (const State *BHD,
                         Vec fg2_fft[3], Vec g, real rho, Vec dg_help);

