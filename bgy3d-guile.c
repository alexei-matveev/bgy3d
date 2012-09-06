/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
#include <libguile.h>
#include <stdbool.h>
#include "petscda.h"            /* PetscInitialize */

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3dH2O.h"
#include "bgy3dH2OS.h"
#include "bgy3d-guile.h"

static char helptext[] = "BGY3d Guile.\n";

/* Update  an  integer  with   the  entry  from  an  association  list
   [intent(in)]  or leave it  unchanged if  there is  no corresponding
   entry. */
static bool alist_getopt_int (SCM alist, const char *key, int *val)
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  bool test = scm_is_true (kv);

  if (test)
    *val = scm_to_int (scm_cdr (kv));

  return test;
}

/* Update  a real  number  with  the entry  from  an association  list
   [intent(in)]  or leave it  unchanged if  there is  no corresponding
   entry: */
static bool alist_getopt_real (SCM alist, const char *key, double *val)
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  bool test = scm_is_true (kv);

  if (test)
    *val = scm_to_double (scm_cdr (kv));

  return test;
}

static ProblemData problem_data (SCM alist)
{
  /* This sets defaults, eventually modified from the command line: */
  ProblemData PD = bgy3d_problem_data ();

  /* Overwrite  defaults with  the  data provided  in the  association
     list. First real (double fprecision) options: */
  alist_getopt_real (alist, "rho", &PD.rho);
  alist_getopt_real (alist, "beta", &PD.beta);
  alist_getopt_real (alist, "norm-tol", &PD.norm_tol);
  alist_getopt_real (alist, "lambda", &PD.lambda);
  alist_getopt_real (alist, "damp-start", &PD.damp);
  alist_getopt_real (alist, "zpad", &PD.zpad);

  real length;
  if (alist_getopt_real (alist, "L", &length)) {
    PD.interval[0] = -length;
    PD.interval[1] = +length;
  }

  /* Integer options: */
  int n;
  if (alist_getopt_int (alist, "N", &n)) {
    PD.N[0] = n;
    PD.N[1] = n;
    PD.N[2] = n;
  }

  alist_getopt_int (alist, "max-iter", &PD.max_iter);
  alist_getopt_int (alist, "solute", &PD.solute);

  /* Preserve invariants or get rid of redundancies: */
  PD.N3 = PD.N[0] * PD.N[1] * PD.N[2];

  length = PD.interval[1] - PD.interval[0];
  PD.h[0] = length / PD.N[0];
  PD.h[1] = length / PD.N[1];
  PD.h[2] = length / PD.N[2];

  return PD;
}

static SCM guile_run_solvent (SCM alist)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  const ProblemData PD = problem_data (alist);

  /* This writes to the disk: */
  BGY3d_solve_2site (&PD, NULL);

  return alist;
}

static SCM guile_run_solute (SCM alist)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  const ProblemData PD = problem_data (alist);

  /* This takes  part of  the input  from the disk,  and writes  to it
     too: */
  BGY3dM_solve_H2O_2site (&PD, NULL);

  return alist;
}

static void inner_main (void *closure, int argc, char **argv)
{
  (void) closure;               /* FIXME: otherwise unused. */

  /*
   * Calling this  will define a  few bgy3d-* gsubrs defined  in these
   * sources:
   */
  scm_c_define_gsubr ("bgy3d-run-solvent", 1, 0, 0, guile_run_solvent);
  scm_c_define_gsubr ("bgy3d-run-solute", 1, 0, 0, guile_run_solute);

  scm_shell (argc, argv);     /* never returns */
}

int bgy3d_guile_main (int argc, char **argv)
{
  /* MPI may choose to rewrite the command line, do it early: */
  PetscInitialize (&argc, &argv, (char*) 0, helptext);

  /* Petsc does not rewrite argv.   Guile will not understand the some
     flags. */
  /* scm_boot_guile (0, NULL, inner_main, NULL); */
  scm_boot_guile (argc, argv, inner_main, NULL);
  return 0; /* never reached */
}
