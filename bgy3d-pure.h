/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void bgy3d_compute_g (Vec g, Vec u0, Vec du);

void bgy3d_pair (State *BHD,
                 const Site a, const Site b, /* struct by value? */
                 Vec f_short[3], Vec f_long[3],
                 Vec u_ini, Vec c2,
                 Vec u2, Vec u2_fft,
                 real damp, real damp_LJ);

void bgy3d_omega (const ProblemData *PD, const DA dc, real rab, Vec w_fft);
void bgy3d_nssa_intra_log (State *BHD, Vec ga_fft, Vec wab_fft, Vec gb, Vec du);

Vec BGY3d_solve_2site (const ProblemData *PD, Vec g_ini);
Vec BGY3d_solve_3site (const ProblemData *PD, Vec g_ini);
void bgy3d_solve_solvent (const ProblemData *PD, int m, const Site solvent[m]);
