/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* Site */
#include "bgy3d-solvents.h"

void bgy3d_sites_dist_mat (int m, const Site sites[m], real r[m][m])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      r[i][j] = sqrt (SQR (sites[i].x[0] - sites[j].x[0]) +
                      SQR (sites[i].x[1] - sites[j].x[1]) +
                      SQR (sites[i].x[2] - sites[j].x[2]));
}

/*
  Prints a table like this:

Solvent:
#       site    ε               σ                q
1       h       0.039710        2.735000         0.200000
2       o       0.514340        3.353000        -0.200000
*/
void bgy3d_sites_show (const char *name, int m, const Site sites[m])
{
  PRINTF ("%s:\n", name);
  PRINTF ("#\tsite\tε       \tσ       \t q\n"); /* unicode here! */
  for (int i = 0; i < m; i++)
    PRINTF ("%d\t%-4s\t%f\t%f\t% f\n",
                 i + 1,
                 sites[i].name,
                 sites[i].epsilon,
                 sites[i].sigma,
                 sites[i].charge);
}

