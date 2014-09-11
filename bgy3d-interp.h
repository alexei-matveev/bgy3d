/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2014 Alexei Matveev
*/

#include "bgy3d.h"              /* DA */

void bgy3d_interp (DA da, const Vec V,      /* real, in */
                   int np, double r[np][3], /* in */
                   double v[np]);           /* out */
