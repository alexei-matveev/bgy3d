/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */

typedef struct Site {
  char name[5];              /* atom types. What are they used for? */
  real x[3];                 /* coordinates */
  real sigma;                /* sigma for LJ */
  real epsilon;              /* epsilon for LJ */
  real charge;               /* charge */
} Site;

/* Get solute sites and name by its name: */
void bgy3d_solute_get (int *n, const Site **sites);

/* Fill intent(out) us[] and uc  fields with the solute field on every
   solvent site. The rest is intent(in). */
void bgy3d_solute_field (const State *BHD,
                         int m, const Site solvent[m], /* m == 2 */
                         Vec us[m], Vec uc, /* intent(out) */
                         int n, const Site solute[n], /* n arbitrary */
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         real damp, real damp_LJ);

