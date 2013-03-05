/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

void hnc3d_solvent_solve (const ProblemData *PD, Vec g[1][1]);
void hnc3d_solute_solve (const ProblemData *PD,
                         int m, const Site solvent[m],
                         int n, const Site solute[n],
                         Vec g[m]);

Vec HNC3d_solvent_solve (const ProblemData *PD, Vec g_ini);
Vec HNC3d_solute_solve (const ProblemData *PD, Vec g_ini);

