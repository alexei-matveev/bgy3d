
/* Details of the struct are left to implementation: */
typedef struct Context Context;

Context* bgy3d_pot_create (const State *BHD, Vec v);
void bgy3d_pot_interp (Context *s, int n, /* const */ real x[n][3], real v[n]);
bool bgy3d_pot_get_value (Context *s, int n, real x[n][3], real v[n], int *p);
void bgy3d_pot_destroy (Context *s);
void bgy3d_pot_test (const State *BHD, Vec vec);
Context* info (const State *BHD,
               int m, const Site solvent[m],
               int n, const Site solute[n],
               Vec g[m],           /* in */
               Vec uc, Vec uc_rho); /* in */

