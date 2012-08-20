real Coulomb_short (real r, real SQRq);
real Coulomb_short_grad (real r, real rx, real SQRq);

real Coulomb_long (real r, real q2);
real Coulomb_long_grad (real r, real rx, real q2);

real Coulomb (real r, real q2);
real Coulomb_grad (real r, real rx, real q2);

void Zeropad_Function (const State *BHD, Vec g, real ZP, real shift);

void ComputeH2O_g (Vec g, Vec g0, Vec dg);

void ImposeBoundaryCondition_Initialize (State *BHD, real zpad);

void Smooth_Function (State *BHD, Vec g, real RL, real RR, real shift);

void ComputeFFTfromCoulomb (State *BHD, Vec uc, Vec f_l[3], fftw_complex *fft_data, real q2, real damp);

void RecomputeInitialData (State *BHD, real damp, real damp_LJ);

void Compute_dg_H2O_inter (State *BHD,
                           Vec f1[3], Vec f1_l[3], Vec g1a, Vec g1b,
                           fftw_complex *coul1_fft, real rho1, real shift1,
                           Vec f2[3], Vec f2_l[3], Vec g2a, Vec g2b,
                           fftw_complex *coul2_fft, real rho2, real shift2,
                           Vec dg, Vec dg_help);

void Compute_dg_H2O_intra (State *BHD, Vec f[3], Vec f_l[3], Vec g1, Vec g2,
                           fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help);

void Compute_dg_H2O_intraIII (State *BHD, Vec f[3], Vec f_l[3], Vec g1, Vec tg,
                             fftw_complex *coul_fft, real rab, Vec dg, Vec dg_help);

void Compute_dg_H2O_intra_ln (State *BHD, Vec g, real rab, Vec dg, Vec dg_help);

void Compute_dg_H2O_normalization_intra (const State *BHD, Vec g, real rab,
                                         Vec dg, Vec dg_help);

void Solve_NormalizationH2O_small (const State *BHD, Vec gc, real rc, Vec g, Vec t,
                                   Vec dg, Vec dg_help, real zpad);
void Solve_NormalizationH2O_smallII (const State *BHD, Vec gc, real rc, Vec g, Vec t,
                                     Vec dg, Vec dg_help, real zpad);

Vec BGY3d_solve_2site (ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3d_solve_3site (ProblemData *PD, Vec g_ini, int vdim);
