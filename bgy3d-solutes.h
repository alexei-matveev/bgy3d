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

/* Fill  fields in the  State struct  with the  solute field  on every
   solvent site: */
void bgy3d_solute_field (State *BHD, int n, const Site S[n], real damp, real damp_LJ);
