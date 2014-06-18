/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2014 Alexei Matveev
*/

/*
  This    one     defines    a    few     scheme    subroutines    and
  returns. Implementation in bgy3d-guile.c.
*/
void bgy3d_guile_init (int argc, char **argv);

/*
  This computes and/or fills the energy *e and gradient g[][] with the
  derivative of solvation term with respect to coordiantes x[][].  See
  bgy3d-guile.c.
*/
void bgy3d_molmech (int n, double x[n][3], double *e, double g[n][3]);
