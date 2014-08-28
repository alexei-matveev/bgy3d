/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

#include <libguile.h>           /* FIXME: SCM */

/* FIXME: Implementation is in bgy3d-guile.c */
bool bgy3d_getopt_test (const char key[]);
bool bgy3d_getopt_bool (const char key[], bool *val);
bool bgy3d_getopt_int (const char key[], int *val);
bool bgy3d_getopt_real (const char key[], double *val);
bool bgy3d_getopt_string (const char key[], int len, char val[len]);
bool bgy3d_getopt_scm (const char key[], SCM *val);
