/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/
#include <libguile.h>

typedef struct Site {
  char name[5];              /* atom types. What are they used for? */
  real x[3];                 /* coordinates */
  real sigma;                /* sigma for LJ */
  real epsilon;              /* epsilon for LJ */
  real charge;               /* charge */
  SCM site; /* Scheme representation. FIXME: should be make it void*? */
} Site;

/* Fill intent(out) us[] and uc  fields with the solute field on every
   solvent site. The rest is intent(in). */
void bgy3d_solute_field (const State *BHD,
                         int m, const Site solvent[m], /* in */
                         int n, const Site solute[n],  /* in */
                         Vec us[m],                    /* out */
                         Vec uc_fft,           /* out, complex */
                         Vec uc_rho,           /* out, optional */
                         Vec uc,               /* out, optional */
                         void (*density)(int k, const real x[k][3], real rho[k]));


/* Differential of solute field with respect to solute coordinates. */
void bgy3d_solute_field1 (const State *BHD,
                          int m, const Site solvent[m], /* in */
                          int n, const Site solute[n],  /* in */
                          real dx[n][3],                /* in */
                          Vec dv[m]);                   /* out */

void bgy3d_solute_form (const State *BHD, int n, const Site solute[n],
                        int m, Vec v_fft[m]); /* inout */


/* Make functions  used as force field primitives  available to Scheme
   code: */
void bgy3d_force_field_init ();
