/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/* Implemented in Fortran. See rism.f90 */
int rism_nrad (const ProblemData *PD);
real rism_rmax (const ProblemData *PD);
ProblemData rism_upscale (const ProblemData *PD);

void rism_solvent (const ProblemData *PD,
                   int m, const Site solvent[m],
                   real t[m][m][*],  /* [m][m][nrad], out */
                   real x[m][m][*]); /* [m][m][nrad], out */

void rism_solute (const ProblemData *PD,
                  int n, const Site solute[n],
                  int m, const Site solvent[m]);
