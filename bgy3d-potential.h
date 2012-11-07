
void bgy3d_pot_alloc (DA da, const ProblemData *PD, Vec v);
void* bgy3d_pot_ini (void);
int bgy3d_pot_get_value (void *s, int n, real x[n][3], real v[n]);
void bgy3d_pot_destroy (void *s);
void bgy3d_pot_test (void);
