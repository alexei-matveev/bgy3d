/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/* Allocates and initializes  a matrix of intra-molecular correlations
   except of diagonal elements that are implicitly 1: */
void bgy3d_omega_fft_create (const State *BHD, int m, const Site solvent[m],
                             Vec omega_fft[m][m]); /* out, creates them */

void bgy3d_nssa_intra_log (State *BHD, Vec ga_fft, Vec wab_fft, Vec gb, Vec du);

void bgy3d_solve_solvent (const ProblemData *PD, int m, const Site solvent[m], Vec g[m][m]);
