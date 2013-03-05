/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

typedef struct Site {
  char name[5];              /* atom types. What are they used for? */
  real x[3];                 /* coordinates */
  real sigma;                /* sigma for LJ */
  real epsilon;              /* epsilon for LJ */
  real charge;               /* charge */
} Site;

/* Get solute sites and name by its name: */
void bgy3d_solute_get (const char *name, int *n, const Site **sites);

/* Fill intent(out) us[] and uc  fields with the solute field on every
   solvent site. The rest is intent(in). */
void bgy3d_solute_field (const State *BHD,
                         int m, const Site solvent[m], /* in */
                         int n, const Site solute[n],  /* in */
                         Vec us[m],                    /* out */
                         Vec uc, Vec uc_rho,           /* out, optional */
                         void (*density)(int k, const real x[k][3], real rho[k]));

