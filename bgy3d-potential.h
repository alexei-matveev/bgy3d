
/* Details of the struct are left to implementation: */
typedef struct Context Context;

Context* bgy3d_pot_create (DA da, const ProblemData *PD, Vec v);
int bgy3d_pot_get_value (Context *s, int n, real x[n][3], real v[n]);
void bgy3d_pot_destroy (Context *s);
void bgy3d_pot_test (Context *s);
