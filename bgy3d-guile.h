/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/
#include <libguile.h>           /* SCM, scm_* */

/* This  one defines  a few  scheme funcitons,  starts the  shell, and
   never returns: */
int bgy3d_guile_main (int argc, char **argv);


static inline void
to_int1 (SCM x, int n, int y[n])
{
  for (int i = 0; i < n; i++)
    {
      y[i] = scm_to_int (scm_car (x));
      x = scm_cdr (x);
    }
}


static inline void
to_double1 (SCM x, int n, double y[n])
{
  for (int i = 0; i < n; i++)
    {
      y[i] = scm_to_double (scm_car (x));
      x = scm_cdr (x);
    }
}


static inline SCM
from_double1 (int n, double x[n])
{
  SCM y = SCM_EOL;

  while (n-- > 0)
    y = scm_cons (scm_from_double (x[n]), y);

  return y;
}


static inline void
to_double2 (SCM x, int m, int n, double y[m][n])
{
  for (int i = 0; i < m; i++)
    {
      to_double1 (scm_car (x), n, y[i]);
      x = scm_cdr (x);
    }
}


static inline SCM
from_double2 (int m, int n, double x[m][n])
{
  SCM y = SCM_EOL;

  while (m-- > 0)
    y = scm_cons (from_double1 (n, x[m]), y);

  return y;
}
