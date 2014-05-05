/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/*==========================================*/
/* Hinweise zur Faltung: */
/* g(x) muss immer um 0 herum bekannt sein. Also sollten auch */
/* die x_M um 0 herum liegen. */
/*===============================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-solutes.h"      /* FIXME: struct Site in prototypes */
#include "bgy3d-pure.h"         /* BGY3d_solvent_solve(), ... */
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-impure.h"       /* BGY3d_solute_solve() */
#include "bgy3d-vec.h"          /* bgy3d_vec_*() */
#include "hnc3d.h"              /* hnc3d_solute_solve(), ... */

#ifdef WITH_GUILE
#include "bgy3d-guile.h"        /* bgy3d_guile_main */
#endif

#ifdef WITH_GUILE
int main (int argc, char **argv)
{
  bgy3d_guile_main (argc, argv);
}
#else
static char helptext[] = "Solving BGY3d equation.\n";

int main (int argc, char **argv)
{
  verbosity = 0;                /* global var */

  /* Petsc  or MPI  may  choose to  rewrite  the command  line, do  it
     early: */
  PetscInitialize (&argc, &argv, (char*)0, helptext);

  /* Make Petsc abort when it encounters an error: */
  PetscPushErrorHandler (PetscAbortErrorHandler, NULL);

  /*
    Read the  command line options.  Petsc insists on keys  having the
    leading dash, so  keep them for the moment.   Set global verbosity
    early  enough. This is  the only  short option!   Use long-options
    prefixed by "--" as the usual convention.
  */
  bgy3d_getopt_int ("-v", &verbosity);

  /* This calls bgy3d_getopt_*() a few more times: */
  ProblemData PD = bgy3d_problem_data ();

  PetscViewerSetFormat(PETSC_VIEWER_STDERR_WORLD, PETSC_VIEWER_ASCII_MATLAB);
  PetscViewerSetFormat(PETSC_VIEWER_STDERR_SELF, PETSC_VIEWER_ASCII_MATLAB);
  /*==================================*/

  {
    int np;
    MPI_Comm_size (PETSC_COMM_WORLD, &np);
    PRINTF ("NP= %d\n", np);
  }

  /* FIXME:  some solvers do  this, some  dont. Sometimes  this output
     appears twice: */
  bgy3d_problem_data_print (&PD);

  /*
    Read method to solve from  command line.  There seem to be several
    solvers that address a  problem of finding solvent distribution in
    external field  given the direct correlation function  of the pure
    solvent. This is a common entry  point for all of them. The actual
    algorith will be affected by --closure HNC/PY/KH and --snes-solver
    newton/picard/jager.

    FIXME:  need   to  parse   database  to  get   the  solvent/solute
    descriptions from  their names.  Currently the  database is parsed
    with   Guile,  so   you   may  want   to   try  WITH_GUILE   build
    instead. Alternatively one  could use Guile for bgy3d_get_solute()
    only. Yet another alternative is  to use easily parseable files as
    inputs.
  */
#error "Needs some work"

  if (bgy3d_getopt_test ("--hnc"))
    {
      if (bgy3d_getopt_test ("--solute"))
        {}                      /* solute/solvent by HNC */
      else
        {}                      /* pure solvent by HNC */
    }

  if (bgy3d_getopt_test ("--bgy"))
    {
      if (bgy3d_getopt_test ("--solute"))
        {}                      /* solvent/solute by BGY */
      else
        {}                      /* pure solvent by BGY */
    }

  int ierr = PetscFinalize();CHKERRQ(ierr);

  /*
   * Original version of BGY3D  executable returned 1. Make interprets
   * non-zero codes as a failure, so this was changed to 0:
   */
  return 0;
}
#endif  /* else of ifdef WITH_GUILE */
