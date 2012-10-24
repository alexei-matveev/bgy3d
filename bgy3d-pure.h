/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
real Coulomb_short (real r, real SQRq);
real Coulomb_short_grad (real r, real rx, real SQRq);

real Coulomb_long (real r, real q2);
real Coulomb_long_grad (real r, real rx, real q2);

real Coulomb (real r, real q2);
real Coulomb_grad (real r, real rx, real q2);

void Zeropad_Function (const State *BHD, Vec g, real ZP, real shift);

void ComputeH2O_g (Vec g, Vec g0, Vec dg);

void Smooth_Function (State *BHD, Vec g, real RL, real RR, real shift);

void ComputeFFTfromCoulomb (State *BHD, Vec uc, Vec f_l[3], Vec uc_fft, real factor);

void RecomputeInitialData (State *BHD, real damp, real damp_LJ);

void Compute_dg_H2O_inter (State *BHD,
                           Vec f1[3], Vec f1_l[3], Vec g1a, Vec g1b,
                           Vec coul1_fft, real rho1,
                           Vec f2[3], Vec f2_l[3], Vec g2a, Vec g2b,
                           Vec coul2_fft, real rho2,
                           Vec dg, Vec dg_help);

#ifdef INTRA1
/* This build seems to be broken: */
void Compute_dg_H2O_intra (State *BHD, Vec f[3], Vec f_l[3], Vec g1, Vec g2,
                           Vec coul_fft, real rab, Vec dg, Vec dg_help);
#else
void Compute_dg_H2O_intraIII (State *BHD, Vec f[3], Vec f_l[3], Vec g1, Vec tg,
                             Vec coul_fft, real rab, Vec dg, Vec dg_help);
#endif

void Compute_dg_H2O_intra_ln (State *BHD, Vec g, real rab, Vec dg);

void Compute_dg_H2O_normalization_intra (const State *BHD, Vec g, real rab,
                                         Vec dg, Vec dg_help);

void Solve_NormalizationH2O_small (const State *BHD, Vec gc, real rc, Vec g, Vec t,
                                   Vec dg, Vec dg_help, real zpad);
void Solve_NormalizationH2O_smallII (const State *BHD, Vec gc, real rc, Vec g, Vec t,
                                     Vec dg, Vec dg_help, real zpad);

void bgy3d_solve_normalization (const State *BHD,
                                Vec gc_fft, /* complex, intent(in) */
                                real rc,
                                Vec g,  /* real, intent(in) */
                                Vec t); /* real, intent(out) */

Vec BGY3d_solve_2site (const ProblemData *PD, Vec g_ini);
Vec BGY3d_solve_3site (const ProblemData *PD, Vec g_ini);
