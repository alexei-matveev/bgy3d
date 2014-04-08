/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/* This  one defines  a few  scheme funcitons,  starts the  shell, and
   never returns: */
int bgy3d_guile_main (int argc, char **argv);

/* This one defines a few scheme subroutines and returns: */
void bgy3d_guile_init (int argc, char **argv);

/* Get solute sites and name by its name: */
void bgy3d_solute_get (const char *name, int *n, Site **sites);
