/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void bgy3d_omega (const ProblemData *PD, const DA dc, real rab, Vec w_fft);
void bgy3d_nssa_intra_log (State *BHD, Vec ga_fft, Vec wab_fft, Vec gb, Vec du);

void bgy3d_solve_solvent (const ProblemData *PD, int m, const Site solvent[m]);

Vec BGY3d_solvent_solve (const ProblemData *PD, Vec g_ini);
Vec BGY3d_solvent_solve_h2o (const ProblemData *PD, Vec g_ini);
