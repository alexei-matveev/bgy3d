/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
#include <libguile.h>
#include <stdbool.h>
#include "petscda.h"            /* PetscInitialize */

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3dH2O.h"
#include "bgy3dH2OS.h"          /* bgy3d_solve_with_solute */
#include "bgy3d-getopt.h"       /* bgy3d_save_vec, bgy3d_load_vec */
#include "bgy3d-guile.h"



#ifndef SIZEOF_VOID_PTR_EQ_4 /* on 32 bit platforms set this CPP var. */
# define scmx_ptr(x) scm_from_uint64 ((intptr_t) (x)) /* scm eXtension */
# define void_ptr(x) ((void*) scm_to_uint64 (x))
#else
# define scmx_ptr(x) scm_from_uint32 ((intptr_t) (x))
# define void_ptr(x) ((void*) scm_to_uint32 (x))
#endif

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
  if (alist_getopt_real (alist, "L", &length))
    {
      PD.interval[0] = -length;
      PD.interval[1] = +length;
    }

  /* Integer options: */
  int n;
  if (alist_getopt_int (alist, "N", &n))
    {
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


/* A Site is represented by an sexp like:

   ("H" (0.6285 0.0 0.0) 2.735 0.03971 0.2) */
static Site make_site (SCM s)
{
  Site S;

  SCM name = scm_car (s);       /* site name */
  SCM x = scm_cadr (s);         /* position */
  SCM ff = scm_cddr (s);        /* force field params */

  /* Currently  max_len  ==  5,  that  is  enough  for  4  letters  ++
     '\0'. Anything that is longer will be truncated: */
  const size_t max_len = sizeof(S.name);

  /* Does not null-terminate, returns the number of bytes required for
     the string not counting the terminating \0: */
  size_t len = scm_to_locale_stringbuf (name, S.name, max_len);

  if (len < max_len)
    S.name[len] = '\0';
  else
    S.name[max_len - 1] = '\0'; /* FIXME: truncation here! */

  // printf ("len=%ld, str=>%s<\n", len, S.name);

  S.x[0] = scm_to_double (scm_car (x));
  S.x[1] = scm_to_double (scm_cadr (x));
  S.x[2] = scm_to_double (scm_caddr (x));

  S.sigma = scm_to_double (scm_car (ff));
  S.epsilon = scm_to_double (scm_cadr (ff));
  S.charge = scm_to_double (scm_caddr (ff));

  return S;
}


static void make_solute (SCM solute, int *n, Site **sites, char **name)
{
  SCM s_name = scm_car (solute);   /* string */
  SCM l_sites = scm_cadr (solute); /* list */
  int length = scm_to_int (scm_length (l_sites));

  /* will you free() it? */
  Site *new = (Site*) malloc (length * sizeof(Site));

  for (int i = 0; i < length; i++)
    {
      new[i] = make_site (scm_car (l_sites));
      l_sites = scm_cdr (l_sites);
    }

  *n = length;
  *sites = new;
  *name = scm_to_locale_string (s_name); /* will you free() it? */
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


static SCM guile_run_solute (SCM solute, SCM settings)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  const ProblemData PD = problem_data (settings);

  int n;
  Site *sites;
  char *name;
  Vec g[2];

  /* Allocates sites, name: */
  make_solute (solute, &n, &sites, &name);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solute is %s.\n", name);

  /* This  takes part  of the  input  from the  disk, returns  solvent
     distribution in Vec g[] (dont forget to destroy them): */
  bgy3d_solve_with_solute (&PD, n, sites, g);

  free (name);
  free (sites);

  /* Caller, dont forget to destroy them! */
  return scm_list_2 (scmx_ptr (g[0]), scmx_ptr (g[1]));
}



static SCM guile_vec_destroy (SCM vec)
{
  PetscErrorCode err = VecDestroy (void_ptr (vec));
  assert (!err);

  return scm_from_int (err);
}


static SCM guile_vec_save (SCM path, SCM vec)
{
  char *c_path = scm_to_locale_string (path); /* free() it! */

  bgy3d_save_vec (c_path, void_ptr (vec));

  free (c_path);

  return SCM_UNDEFINED;
}


static SCM guile_vec_load (SCM path)
{
  char *c_path = scm_to_locale_string (path); /* free() it! */

  Vec vec = bgy3d_load_vec (c_path);

  free (c_path);

  return scmx_ptr (vec);
}


static SCM guile_vec_length (SCM vec)
{
  int len;
  Vec c_vec = void_ptr (vec);

  VecGetSize (c_vec, &len);

  return scm_from_int (len);
}


static SCM guile_vec_ref (SCM vec, SCM ix)
{
  real vals[1];

  Vec c_vec = void_ptr (vec);

  int keys[1] = {scm_to_int (ix)};

  VecGetValues (c_vec, 1, keys, vals);

  return scm_from_double (vals[0]);
}


static SCM guile_test (SCM inp)
{
  Vec v = void_ptr (inp);

  return scmx_ptr (v);
}


/* Finalize Petsc and MPI */
static void finalize (void)
{
  PetscErrorCode err = PetscFinalize ();
  assert (!err);
}



void bgy3d_guile_init (int argc, char **argv)
{
  assert (sizeof (Vec) == sizeof (void*));
#ifndef SIZEOF_VOID_PTR_EQ_4
  assert (sizeof (void*) == sizeof (uint64_t));
#else
  assert (sizeof (void*) == sizeof (uint32_t));
#endif

  /* MPI may  choose to rewrite the  command line, do  it early. Petsc
     does not rewrite argv.  Guile will not understand Petsc flags. */
  PetscInitialize (&argc, &argv, NULL, helptext);

  /* Add  an  exit handler  that  calls  PetscFinalize(). Executed  by
     exit() according to POSIX: */
  atexit (finalize);

  /*
   * Calling this  will define a  few bgy3d-* gsubrs defined  in these
   * sources:
   */
  scm_c_define_gsubr ("bgy3d-run-solvent", 1, 0, 0, guile_run_solvent);
  scm_c_define_gsubr ("bgy3d-run-solute", 2, 0, 0, guile_run_solute);
  scm_c_define_gsubr ("bgy3d-vec-destroy", 1, 0, 0, guile_vec_destroy);
  scm_c_define_gsubr ("bgy3d-vec-save", 2, 0, 0, guile_vec_save);
  scm_c_define_gsubr ("bgy3d-vec-load", 1, 0, 0, guile_vec_load);
  scm_c_define_gsubr ("bgy3d-vec-length", 1, 0, 0, guile_vec_length);
  scm_c_define_gsubr ("bgy3d-vec-ref", 2, 0, 0, guile_vec_ref);
  scm_c_define_gsubr ("bgy3d-test", 1, 0, 0, guile_test);
}

static void inner_main (void *closure, int argc, char **argv)
{
  (void) closure;               /* FIXME: otherwise unused. */

  bgy3d_guile_init (argc, argv);

  scm_shell (argc, argv);     /* never returns */
}

int bgy3d_guile_main (int argc, char **argv)
{
  /* Petsc does not rewrite argv.   Guile will not understand the some
     flags. */
  /* scm_boot_guile (0, NULL, inner_main, NULL); */
  scm_boot_guile (argc, argv, inner_main, NULL);
  return 0; /* never reached */
}
