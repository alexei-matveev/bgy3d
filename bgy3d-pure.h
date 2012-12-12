/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
real Coulomb_short (real r, real SQRq);
real Coulomb_short_grad (real r, real rx, real SQRq);

real Coulomb_long (real r, real q2);
real Coulomb_long_grad (real r, real rx, real q2);

real Coulomb (real r, real q2);
real Coulomb_grad (real r, real rx, real q2);

void ComputeH2O_g (Vec g, Vec g0, Vec dg);

void ComputeFFTfromCoulomb (State *BHD,
                            Vec uc, Vec fc[3], /* intent(out) */
                            Vec uc_fft,    /* complex, intent(out) */
                            Vec fc_fft[3], /* complex, intent(out) */
                            real factor);


void Compute_dg_H2O_intra_ln (State *BHD, Vec g, real rab, Vec dg);

void bgy3d_solve_normalization (const State *BHD,
                                Vec gc_fft, /* complex, intent(in) */
                                real rc,
                                Vec g,  /* real, intent(in) */
                                Vec t); /* real, intent(out) */

Vec BGY3d_solve_2site (const ProblemData *PD, Vec g_ini);
Vec BGY3d_solve_3site (const ProblemData *PD, Vec g_ini);

/*
  Does the mixing:

    dg := a * dg_new + (1 - a) * dg

  Returns the norm of the difference |dg_new - dg|.
 */
real bgy3d_vec_mix (Vec dg, Vec dg_new, real a, Vec work);
