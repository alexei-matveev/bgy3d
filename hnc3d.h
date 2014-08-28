/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/
#include <libguile.h>           /* SCM */

void hnc3d_solvent_solve (const ProblemData *PD,
                          int m, const Site solvent[m],
                          Vec g[m][m]);

void hnc3d_solute_solve (const ProblemData *PD,
                         const int m, const Site solvent[m],
                         const int n, const Site solute[n],
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         SCM *dict, /* inout, association list */
                         Vec g[m],  /* out */
                         Context **medium,   /* out */
                         Restart **restart); /* inout */
