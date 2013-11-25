/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include <libguile.h>           /* SCM, SCM_EOL */


/* Empty association  list will be  accepted as the first  argument of
   rism_solvent() and rism_solute(): */
#define RISM_NULL_ENV SCM_EOL

/* Implemented in Fortran. See rism.f90 */
int rism_nrad (const ProblemData *PD);
real rism_rmax (const ProblemData *PD);
ProblemData rism_upscale (const ProblemData *PD);
double rism_self_energy (int n, const Site sites[n], const int spec[n]);

void rism_solvent (SCM env,     /* SCM alist */
                   const ProblemData *PD,
                   int m, const Site solvent[m],
                   real t[m][m][*],  /* [m][m][nrad] or NULL, out */
                   real x[m][m][*],  /* [m][m][nrad] or NULL, out */
                   void *retval);    /* SCM* or NULL, out */

void rism_solute (SCM env,      /* SCM alist */
                  const ProblemData *PD,
                  int n, const Site solute[n],
                  int m, const Site solvent[m],
                  real x[m][m][*], /* [m][m][nrad] or NULL, in */
                  void *retval);   /* SCM* or NULL, out */
