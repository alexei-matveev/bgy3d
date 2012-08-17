int bgy3d_getopt_test (const char key[]);
int bgy3d_getopt_int (const char key[], int *val);
int bgy3d_getopt_real (const char key[], double *val);

void bgy3d_load_vec (const char file[], Vec *vec);
void bgy3d_save_vec (const char file[], const Vec vec);
