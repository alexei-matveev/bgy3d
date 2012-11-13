/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
#include <libguile.h>
#include <stdbool.h>
#include "petscda.h"            /* PetscInitialize */

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-pure.h"
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-impure.h"       /* bgy3d_solve_with_solute */
#include "bgy3d-getopt.h"       /* bgy3d_save_vec, bgy3d_load_vec */
#include "bgy3d-fft.h"          /* bgy3d_fft_test() */
#include "bgy3d-guile.h"



/* We  need  to apply  either  scm_from_uint32() or  scm_from_uint64()
   depending on the size of intptr_t.  This would be a good use of C11
   type generic macros.  SCMX stays  for SCM eXtension in order not to
   confuse the  macros with those in  libguile.h. Here we  make use of
   the macro defined by libguile.h for the size of intptr_t: */
#if SCM_SIZEOF_INTPTR_T == 8
# define scmx_from_intptr scm_from_uint64
# define scmx_to_intptr scm_to_uint64
#elif SCM_SIZEOF_INTPTR_T == 4
# define scmx_from_intptr scm_from_uint32
# define scmx_to_intptr scm_to_uint32
#else
# error "unknown pointer size!"
#endif

#define scmx_ptr(x) scmx_from_intptr ((intptr_t) (x))
#define void_ptr(x) ((void*) scmx_to_intptr (x))

static char helptext[] = "BGY3d Guile.\n";


/* Update an integer with the  entry from an association list or leave
   it  unchanged if  there is  no  corresponding entry.  Alist is  not
   modified. */
static bool alist_getopt_int (SCM alist, const char *key, int *val)
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  bool test = scm_is_true (kv);

  if (test)
    *val = scm_to_int (scm_cdr (kv));

  return test;
}


/* Update a  real number  with the entry  from an association  list or
   leave it unchanged if there is no corresponding entry. Alist is not
   modified. */
static bool alist_getopt_real (SCM alist, const char *key, double *val)
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  bool test = scm_is_true (kv);

  if (test)
    *val = scm_to_double (scm_cdr (kv));

  return test;
}

/* Update a function  pointer with the entry from  an association list
   or leave it unchanged if there is no corresponding entry: */
static bool alist_getopt_funptr (SCM alist, /* intent(in) */
                                 const char *key,
                                 void (**val)(void))
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  bool test = scm_is_true (kv);

  /* FIXME: potentially non-portable conversion of an integer to void*
     and then to function pointer, void (*)(): */
  if (test)
    *val = void_ptr (scm_cdr (kv));

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

  /* This declares and  sets a function pointer. If  the settings dont
     specify it, it should remain NULL: */
  void (*qm_density) (int n, const real x[n][3], real rho[n]) = NULL;

  /* Cast is to silence the warning here.  Note that we pass a pointer
     to  a  funptr,  void  (**)(),  as the  function  is  supposed  to
     (eventually) set that funptr to something meaningful: */
  alist_getopt_funptr (settings, "qm-density", (void (**)()) &qm_density);

  // printf ("qm-density=%p\n", qm_density); /* print funptr in hex */

  /*
    This  takes part  of  the  input from  the  disk, returns  solvent
    distribution  in Vec  g[] (dont  forget to  destroy them).   If no
    additional charge distribution is  associated with the solute pass
    NULL as  the function  pointer. Similarly, if  you do not  want an
    iterator over the solvent potential pass NULL:
  */
  Context *iter;
  bgy3d_solve_with_solute (&PD, n, sites, qm_density, g, &iter);

  free (name);
  free (sites);

  /* Caller, dont forget to destroy them! */
  SCM gs = scm_list_2 (scmx_ptr (g[0]), scmx_ptr (g[1]));
  SCM v = scmx_ptr (iter);

  /* Return multiple values: */
  return scm_values (scm_list_2 (gs, v));
}



static SCM guile_pot_destroy (SCM iter)
{
  bgy3d_pot_destroy (void_ptr (iter));

  return scm_from_int (0);
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

/* Return MPI runk in PETSC_COMM_WORLD: */
static SCM guile_rank (void)
{
  int rank;
  MPI_Comm_rank (PETSC_COMM_WORLD, &rank);
  return scm_from_int (rank);
}


/* Return MPI size of PETSC_COMM_WORLD: */
static SCM guile_size (void)
{
  int size;
  MPI_Comm_size (PETSC_COMM_WORLD, &size);
  return scm_from_int (size);
}


/* Reduce buffer by summing respective entries on all workeres: */
static void comm_allreduce (void *buf, int count, MPI_Datatype type)
{
  int err = MPI_Allreduce (MPI_IN_PLACE, buf, count, type, MPI_SUM,
                           PETSC_COMM_WORLD);
  assert (!err);
}

/* An inefficient way of getting just one value, vec[ix], even if that
   value is not stored locally. Collective. */
static SCM guile_vec_ref (SCM vec, SCM ix)
{
  assert (sizeof(real) == sizeof(double)); /* See MPI_DOUBLE */

  Vec c_vec = void_ptr (vec);
  int i = scm_to_int (ix);

  PetscInt lo, hi;
  VecGetOwnershipRange (c_vec, &lo, &hi);

  PetscInt keys[1] = {i};
  real vals[1];
  if (lo <= i && i < hi)
    VecGetValues (c_vec, 1, keys, vals); /* local */
  else
    vals[0] = 0.0;

  /* Make result known on all workers: */
  comm_allreduce (vals, 1, MPI_DOUBLE);

  return scm_from_double (vals[0]);
}


static SCM guile_test (SCM m, SCM n, SCM k)
{
  return scm_from_double (bgy3d_fft_test (scm_to_int (m),
                                          scm_to_int (n),
                                          scm_to_int (k)));
}


/* Finalize Petsc and MPI */
static void finalize (void)
{
  PetscErrorCode err = PetscFinalize ();
  assert (!err);
}



/* Calling this will define a few bgy3d-* gsubrs introduced above: */
static SCM guile_bgy3d_module_init (void)
{
  /* If Scheme executes this code  inside a module (which we do), then
     all these gsubrs will be  module procedures available only in the
     module itself or by an explicit (use-modules ...): */
  scm_c_define_gsubr ("bgy3d-run-solvent", 1, 0, 0, guile_run_solvent);
  scm_c_define_gsubr ("bgy3d-run-solute", 2, 0, 0, guile_run_solute);
  scm_c_define_gsubr ("bgy3d-pot-destroy", 1, 0, 0, guile_pot_destroy);
  scm_c_define_gsubr ("bgy3d-vec-destroy", 1, 0, 0, guile_vec_destroy);
  scm_c_define_gsubr ("bgy3d-vec-save", 2, 0, 0, guile_vec_save);
  scm_c_define_gsubr ("bgy3d-vec-load", 1, 0, 0, guile_vec_load);
  scm_c_define_gsubr ("bgy3d-vec-length", 1, 0, 0, guile_vec_length);
  scm_c_define_gsubr ("bgy3d-vec-ref", 2, 0, 0, guile_vec_ref);
  scm_c_define_gsubr ("bgy3d-rank", 0, 0, 0, guile_rank);
  scm_c_define_gsubr ("bgy3d-size", 0, 0, 0, guile_size);
  scm_c_define_gsubr ("bgy3d-test", 3, 0, 0, guile_test);

  return SCM_UNSPECIFIED;
}

void bgy3d_guile_init (int argc, char **argv)
{
  /* The  code  above  assumes an  opaque  type  Vec  can be  cast  to
     void*: */
  assert (sizeof (Vec) == sizeof (void*));

  /* Depending on this preprocessor flag,  void* is cast either to 32-
     or 64-bit integers.  Make sure we  use the right size. Would be a
     good use for a static (compile time) assert: */
#if SCM_SIZEOF_INTPTR_T == 8
  assert (sizeof (void*) == sizeof (uint64_t));
#else
  assert (sizeof (void*) == sizeof (uint32_t));
#endif

  /* We assume that  data- and function pointers are  of the same size
     when assigning a void* to void (*)() in alist_getopt_funptr(): */
  {
    void (*fn)(void);           /* only for an assert */
    assert (sizeof (fn) == sizeof (void*));
  }

  /* MPI may  choose to rewrite the  command line, do  it early. Petsc
     does not rewrite argv.  Guile will not understand Petsc flags. */
  PetscInitialize (&argc, &argv, NULL, helptext);

  /* Add  an  exit handler  that  calls  PetscFinalize(). Executed  by
     exit() according to POSIX: */
  atexit (finalize);

  /* Make Petsc abort when it encounters an error: */
  PetscPushErrorHandler (PetscAbortErrorHandler, NULL);

  /*
   * Note  that  the names  defined  here  are  put into  the  private
   * namespace  of (guile-user)  module. If  you want  to call  any of
   * these you may need to "steal" it from there by dereferencing them
   * as e.g. (@@ (guile-user) guile-bgy3d-module-init).
   *
   * Calling this will define bgy3d-* gsubrs:
   */
  scm_c_define_gsubr ("guile-bgy3d-module-init", 0, 0, 0, guile_bgy3d_module_init);
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
