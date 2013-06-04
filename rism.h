/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/* Implemented in Fortran. See rism.f90 */
void rism_solvent (const ProblemData *PD,
                   int m, const Site solvent[m],
                   real chi_fft[m][m][*]); /* [m][m][nrad], out */

void rism_solute (const ProblemData *PD,
                  int n, const Site solute[n],
                  int m, const Site solvent[m]);
