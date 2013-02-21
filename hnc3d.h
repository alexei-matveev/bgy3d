/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

Vec hnc3d_solve_newton (const ProblemData *PD, Vec g_ini);
Vec hnc3d_solve_picard (const ProblemData *PD, Vec g_ini);

Vec hnc3d_solute_solve_picard (const ProblemData *PD, Vec g_ini);
Vec hnc3d_solute_solve_newton (const ProblemData *PD, Vec g_ini);

