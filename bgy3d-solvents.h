/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/* smoothing parameter */
#define G_COULOMB_INVERSE_RANGE 1.2 //1.2

void bgy3d_sites_dist_mat (int m, const Site sites[m], real r[m][m]);
void bgy3d_sites_show (const char *name, int m, const Site sites[m]);
