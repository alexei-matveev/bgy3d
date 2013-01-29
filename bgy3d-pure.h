/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void ComputeH2O_g (Vec g, Vec u0, Vec du);

void ComputeFFTfromCoulomb (State *BHD,
                            Vec uc, Vec fc[3], /* intent(out) */
                            Vec uc_fft,    /* complex, intent(out) */
                            Vec fc_fft[3], /* complex, intent(out) */
                            real factor);

void bgy3d_omega (const ProblemData *PD, const DA dc, real rab, Vec w_fft);
void bgy3d_nssa_intra_log (State *BHD, Vec ga_fft, Vec wab_fft, Vec gb, Vec du);

Vec BGY3d_solve_2site (const ProblemData *PD, Vec g_ini);
Vec BGY3d_solve_3site (const ProblemData *PD, Vec g_ini);
