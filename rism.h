/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013, 2014 Alexei Matveev
*/

/* Implemented in Fortran. See rism.f90 */
int rism_nrad (const ProblemData *PD);
real rism_rmax (const ProblemData *PD);
ProblemData rism_upscale (const ProblemData *PD);

SCM rism_self_energy (int n, const Site sites[n], const int spec[n]);

void rism_solvent (const ProblemData *PD,
                   int m, const Site solvent[m],
                   real t[m][m][*],  /* [m][m][nrad] or NULL, out */
                   real x[m][m][*],  /* [m][m][nrad] or NULL, out */
                   void *retval);    /* SCM* or NULL, out */

void rism_solute (const ProblemData *PD,
                  int n, const Site solute[n],
                  int m, const Site solvent[m],
                  real x[m][m][*], /* [m][m][nrad] or NULL, in */
                  void *retval);   /* SCM* or NULL, out */

/*
  subroutine rism_solute_renorm &
    (m, solvent, rmax, nrad, x_kvv, alpha, s_kv) bind (c)
*/
void rism_solute_renorm (int m, const Site solvent[m],
                         real rmax, int nrad, real x_fft[m][m][nrad], real G,
                         real s_fft[m][nrad]); /* out */


/*
  Implemented in Fortran. See ./closures.f90:

  subroutine rism_closure (method, beta, n, v, t, c) bind (c)
*/
void rism_closure (int method, real beta,
                   int n, const real v[n], const real t[n], /* in */
                   real c[n]);                              /* out */

