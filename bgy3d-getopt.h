/* -*- mode: c; c-basic-offset: 2; -*- */
int bgy3d_getopt_test (const char key[]);
int bgy3d_getopt_int (const char key[], int *val);
int bgy3d_getopt_real (const char key[], double *val);
int bgy3d_getopt_string (const char key[], char *val, size_t len);

void bgy3d_load_vec (const char file[], Vec *vec);
void bgy3d_save_vec (const char file[], const Vec vec);
void bgy3d_save_vec_ascii (const char file[], const Vec vec);
