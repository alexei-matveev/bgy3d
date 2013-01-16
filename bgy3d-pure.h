/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
real Coulomb_short (real r, real q2);
real Coulomb_short_grad (real r, real rx, real q2);

real Coulomb_long (real r, real q2);
real Coulomb_long_grad (real r, real rx, real q2);

real Coulomb (real r, real q2);
real Coulomb_grad (real r, real rx, real q2);

void ComputeH2O_g (Vec g, Vec u0, Vec du);

void ComputeFFTfromCoulomb (State *BHD,
                            Vec uc, Vec fc[3], /* intent(out) */
                            Vec uc_fft,    /* complex, intent(out) */
                            Vec fc_fft[3], /* complex, intent(out) */
                            real factor);


void Compute_dg_H2O_intra_ln (State *BHD, Vec g, real rab, Vec dg);

void bgy3d_nssa_gamma_cond (const State *BHD, Vec gac_fft, real rbc, Vec gab,
                            Vec tab); /* intent(out) */

Vec BGY3d_solve_2site (const ProblemData *PD, Vec g_ini);
Vec BGY3d_solve_3site (const ProblemData *PD, Vec g_ini);
