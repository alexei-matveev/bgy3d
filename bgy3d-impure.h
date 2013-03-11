/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/*
  A Restart* token to be  passed there an back for resuming iterations
  without  disk access.  There is  no explicit  constructor  for these
  objects     but    they     are    (optionally)     returned    from
  bgy3d_solute_solve().
*/
typedef struct Restart Restart;

/*
  This function is the main entry  point for the BGY3dM equation for a
  m-site solvent and an arbitrary solute.  The vectors in

  Vec g[m], intent(out)

  are  initialzed as  global distributed  arrays and  filled  with the
  solvent site  distributions. It is the responsibility  of the caller
  to destroy them when no more needed.
*/
void bgy3d_solute_solve (const ProblemData *PD,
                         int m, const Site solvent[m],
                         int n, const Site solute[n],
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         Vec g[m],          /* out */
                         Context **medium,  /* out, optional */
                         Restart **restart); /* inout, optional */

void bgy3d_restart_destroy (Restart *restart);

Vec BGY3dM_solve_H2O_2site(const ProblemData *PD, Vec g_ini);

#ifdef WITH_EXTRA_SOLVERS
Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini);
void RecomputeInitialFFTs (State *BHD,
                           int m,
                           Vec g2[m][m],        /* real, in */
                           Vec fs_g2_fft[m][m][3], /* complex, out */
                           Vec fl_g2_fft[m][m][3], /* complex, out */
                           Vec u2[m][m],     /* real, out */
                           Vec u2_fft[m][m]) /* complex, out */
  __attribute__ ((deprecated));

void Compute_H2O_interS (const State *BHD,
                         Vec fg2_fft[3], Vec g, real rho, Vec dg_help)
  __attribute__ ((deprecated));
#endif

