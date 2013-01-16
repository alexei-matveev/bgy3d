/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

int bgy3d_getopt_test (const char key[]);
int bgy3d_getopt_int (const char key[], int *val);
int bgy3d_getopt_real (const char key[], double *val);
int bgy3d_getopt_string (const char key[], char *val, size_t len);
