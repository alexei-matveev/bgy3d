/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */

typedef struct Site {
  char name[5];              /* atom types. What are they used for? */
  real x[3];                 /* coordinates */
  real sigma;                /* sigma for LJ */
  real epsilon;              /* epsilon for LJ */
  real charge;               /* charge */
} Site;

/* Get solute sites and name by index: */
void bgy3d_solute_get (int solute, int *n, const Site **sites, const char **name);

/* Fill  intent(out) us[]  and ul[]  fields with  the solute  field on
   every solvent site. The rest is intent(in). */
void bgy3d_solute_field (const State *BHD, Vec us[2], Vec ul[2],
                         int n, const Site S[n], real damp, real damp_LJ);
