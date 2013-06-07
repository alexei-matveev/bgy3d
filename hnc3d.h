/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void hnc3d_solvent_solve (const ProblemData *PD,
                          int m, const Site solvent[m],
                          Vec g[m][m]);

void hnc3d_solute_solve (const ProblemData *PD,
                         const int m, const Site solvent[m],
                         const int n, const Site solute[n],
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         Vec g[m],
                         Context **medium,   /* out */
                         Restart **restart); /* inout */

Vec HNC3d_solvent_solve (const ProblemData *PD, Vec g_ini);
Vec HNC3d_solute_solve (const ProblemData *PD, Vec g_ini);

