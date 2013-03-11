/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
#include <libguile.h>
#include <stdbool.h>
#include "petscda.h"            /* PetscInitialize */

#include "bgy3d.h"
#include "bgy3d-getopt.h"       /* bgy3d_getopt_test() */
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* bgy3d_solvent_get() */
#include "bgy3d-pure.h"
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-impure.h"       /* bgy3d_solve_with_solute */
#include "hnc3d.h"              /* hnc3d_solute_solve() */
#include "bgy3d-vec.h"          /* bgy3d_vec_save, bgy3d_vec_load */
#include "bgy3d-fft.h"          /* bgy3d_fft_test() */
#include "bgy3d-fftw.h"         /* bgy3d_fft_interp() */
#include "bgy3d-guile.h"



static char helptext[] = "BGY3d Guile.\n";

/*
  We  need  to  apply  either scm_from_uint32()  or  scm_from_uint64()
  depending on the size of intptr_t.   This would be a good use of C11
  type  generic macros.   Here we  make use  of the  macro  defined by
  libguile.h for the size of intptr_t:
*/
#if SCM_SIZEOF_INTPTR_T == 8
# define from_intptr scm_from_uint64
# define to_intptr scm_to_uint64
#elif SCM_SIZEOF_INTPTR_T == 4
# define from_intptr scm_from_uint32
# define to_intptr scm_to_uint32
#else
# error "unknown pointer size!"
#endif

/* The  const qualifier  might be  misleading  with all  the casts  to
   integers and back: */
static SCM from_pointer (const void *ptr)
{
  return from_intptr ((intptr_t) ptr);
}

static void* to_pointer (SCM ptr)
{
  return (void*) to_intptr (ptr);
}

/* Update SCM value  with the entry from an  association list or leave
   it  unchanged if  there is  no corresponding  entry.  Alist  is not
   modified. */
static bool alist_getopt_scm (SCM alist, const char *key, SCM *val)
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  bool test = scm_is_true (kv);

  if (test)
    *val = scm_cdr (kv);

  return test;
}

/* Update an integer with the  entry from an association list or leave
   it  unchanged if  there is  no  corresponding entry.  Alist is  not
   modified. */
static bool alist_getopt_int (SCM alist, const char *key, int *val)
{
  SCM scm;
  bool test = alist_getopt_scm (alist, key, &scm);

  if (test)
    *val = scm_to_int (scm);

  return test;
}


/* Update a  real number  with the entry  from an association  list or
   leave it unchanged if there is no corresponding entry. Alist is not
   modified. */
static bool alist_getopt_real (SCM alist, const char *key, double *val)
{
  SCM scm;
  bool test = alist_getopt_scm (alist, key, &scm);

  if (test)
    *val = scm_to_double (scm);

  return test;
}

/* Update a function  pointer with the entry from  an association list
   or leave it unchanged if there is no corresponding entry: */
static bool alist_getopt_funptr (SCM alist, /* intent(in) */
                                 const char *key,
                                 void (**val)(void))
{
  SCM scm;
  bool test = alist_getopt_scm (alist, key, &scm);

  /* FIXME: potentially non-portable conversion of an integer to void*
     and then to function pointer, void (*)(): */
  if (test)
    *val = to_pointer (scm);

  return test;
}


static ProblemData problem_data (SCM alist)
{
  /* This sets defaults, eventually modified from the command line: */
  ProblemData PD = bgy3d_problem_data ();

  /* Overwrite  defaults with  the  data provided  in the  association
     list. First real (double precision) options: */
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

static void to_array (SCM x, int n, double y[n])
{
  for (int i = 0; i < n; i++)
    {
      y[i] = scm_to_double (scm_car (x));
      x = scm_cdr (x);
    }
}

/* A Site is represented by an sexp like:

   ("H" (0.6285 0.0 0.0) 2.735 0.03971 0.2) */
static Site to_site (SCM s)
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

  /* Fill an array of length three from the list: */
  to_array (x, 3, S.x);

  S.sigma = scm_to_double (scm_car (ff));
  S.epsilon = scm_to_double (scm_cadr (ff));
  S.charge = scm_to_double (scm_caddr (ff));

  return S;
}


static void to_sites (SCM molecule, int *n, Site **sites, char **name)
{
  SCM molecule_name = scm_car (molecule); /* string */
  SCM site_list = scm_cadr (molecule);    /* list */
  const int length = scm_to_int (scm_length (site_list));

  /* will you free() it? */
  Site *new = (Site*) malloc (length * sizeof(Site));

  for (int i = 0; i < length; i++)
    {
      new[i] = to_site (scm_car (site_list));
      site_list = scm_cdr (site_list);
    }

  *n = length;
  *sites = new;                                 /* free() it! */
  *name = scm_to_locale_string (molecule_name); /* free() it! */
}


/* The following  code declares a  State SMOB primarily to  make array
   descriptors, FFT, and laplace matrices available to Scheme: */
static scm_t_bits state_tag;
static scm_t_bits vec_tag;


static SCM from_state (const State *BHD)
{
  SCM obj;
  SCM_NEWSMOB (obj, state_tag, BHD);
  return obj;
}


static SCM from_vec (Vec vec)
{
  SCM obj;
  SCM_NEWSMOB (obj, vec_tag, vec);
  return obj;
}


/* Build a list starting from the tail: */
static SCM from_vec1 (int m, const Vec g[m])
{
  SCM gs = SCM_EOL;             /* empty list */
  for (int i = m - 1; i >= 0; i--)
    gs = scm_cons (from_vec (g[i]), gs);
  return gs;
}


static State* to_state (SCM state)
{
  scm_assert_smob_type (state_tag, state);

  return (State*) SCM_SMOB_DATA (state);
}


static Vec to_vec (SCM vec)
{
  scm_assert_smob_type (vec_tag, vec);

  return (Vec) SCM_SMOB_DATA (vec);
}


static SCM state_make (SCM alist)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  ProblemData *PD = malloc (sizeof *PD); /* free() in state_free() */
  *PD = problem_data (alist);

  return from_state (bgy3d_state_make (PD));
}


static SCM vec_make (SCM state)
{
  State *BHD = to_state (state);
  Vec vec = bgy3d_vec_create (BHD->da);
  return from_vec (vec);
}


static SCM vec_make_complex (SCM state)
{
  State *BHD = to_state (state);
  Vec vec = bgy3d_vec_create (BHD->dc);
  return from_vec (vec);
}


static size_t state_free (SCM state)
{
  State *BHD = to_state (state);
  /*
    If state-destroy was called explicitly  on this smob then the real
    data  is already  freed  and the  pointer  must have  been set  to
    NULL:
  */
  if (BHD)
    {
      free ((void *) BHD->PD);  /* FIXME: discards const! */
      bgy3d_state_destroy (BHD);
    }
  return 0;
}


static size_t vec_free (SCM vec)
{
  Vec c_vec = to_vec (vec);
  if (c_vec)
    bgy3d_vec_destroy (&c_vec);
  return 0;
}


static SCM state_destroy (SCM state)
{
  assert (to_state (state) != NULL);
  state_free (state);
  SCM_SET_SMOB_DATA (state, NULL);
  return SCM_UNSPECIFIED;
}


static SCM vec_destroy (SCM vec)
{
  assert (to_vec (vec) != NULL);
  vec_free (vec);
  SCM_SET_SMOB_DATA (vec, NULL);
  return SCM_UNSPECIFIED;
}


static SCM noop_mark (SCM smob)
{
  /* We dont store SCM values in smobs yet: */
  (void) smob;
  return SCM_BOOL_F;
}


static SCM vec_save (SCM path, SCM vec)
{
  char *c_path = scm_to_locale_string (path); /* free() it! */

  bgy3d_vec_save (c_path, to_vec (vec));

  free (c_path);

  return SCM_UNDEFINED;
}


static SCM vec_load (SCM path)
{
  char *c_path = scm_to_locale_string (path); /* free() it! */

  Vec vec = bgy3d_vec_load (c_path);

  free (c_path);

  return from_vec (vec);
}


static SCM vec_length (SCM vec)
{
  int len;
  VecGetSize (to_vec (vec), &len);

  return scm_from_int (len);
}

/* An inefficient way of getting just one value, vec[ix], even if that
   value is not stored locally. Collective. */
static SCM vec_ref (SCM vec, SCM ix)
{
  assert (sizeof(real) == sizeof(double)); /* See MPI_DOUBLE */

  Vec c_vec = to_vec (vec);
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
  bgy3d_comm_allreduce (1, vals);

  return scm_from_double (vals[0]);
}


static SCM vec_dot (SCM x, SCM y)
{
  return scm_from_double (bgy3d_vec_dot (to_vec (x), to_vec (y)));
}


static SCM vec_fft (SCM state, SCM x)
{
  State *BHD = to_state (state);
  SCM y = vec_make_complex (state);
  Vec y_ = to_vec (y);

  MatMult (BHD->fft_mat, to_vec (x), y_);
  VecScale (y_, 1.0 / sqrt (BHD->PD->N3));

  return y;
}


static SCM vec_ifft (SCM state, SCM y)
{
  State *BHD = to_state (state);
  SCM x = vec_make (state);

  Vec x_ = to_vec (x);
  MatMultTranspose (BHD->fft_mat, to_vec (y), x_);
  VecScale (x_, 1.0 / sqrt (BHD->PD->N3));

  return x;
}


static SCM vec_fft_interp (SCM state, SCM Y, SCM x)
{
  State *BHD = to_state (state);
  Vec Y_ = to_vec (Y);

  double x_[1][3], y[1];
  for (int i = 0; i < 3; i ++)
    {
      x_[0][i] = scm_to_double (scm_car (x));
      x = scm_cdr (x);
    }

  bgy3d_fft_interp (BHD->fft_mat, Y_, 1, x_, y);

  /* FIXME: FFT in BGY3d code is unnormalized: */
  return scm_from_double (y[0] * sqrt (BHD->PD->N3));
}


static SCM vec_set_random (SCM x)
{
  VecSetRandom (to_vec (x), NULL);

  return x;
}

/*
  NOTE:  scm_from_double()   allocates  on  heap   producing  lots  of
  garbage. A  good thing  about it is  that it  also causes GC  to run
  regularly and collect Vec smobs.
*/
static SCM vec_map1 (SCM f, SCM x)
{
  Vec x_ = to_vec (x);
  Vec y_ = bgy3d_vec_duplicate (x_);

  double f_ (double x)
  {
    return scm_to_double (scm_call_1 (f, scm_from_double (x)));
  }
  bgy3d_vec_map1 (y_, f_, x_);

  return from_vec (y_);
}


static SCM vec_map2 (SCM f, SCM x, SCM y)
{
  Vec x_ = to_vec (x);
  Vec y_ = to_vec (y);
  Vec z_ = bgy3d_vec_duplicate (x_);

  double f_ (double x, double y)
  {
    return scm_to_double (scm_call_2 (f, scm_from_double (x), scm_from_double (y)));
  }
  bgy3d_vec_map2 (z_, f_, x_, y_);

  return from_vec (z_);
}


static int state_print (SCM state, SCM port, scm_print_state *pstate)
{
  (void) pstate;
  const State *BHD = to_state (state);
  const ProblemData *PD = BHD->PD;

  scm_puts ("#<state addr: ", port);
  scm_display (from_pointer (BHD), port);
  scm_puts (", shape: ", port);
  for (int i = 0; i < 3; i++)
    {
      scm_display (scm_from_int (PD->N[i]), port);
      if (i != 2)
        scm_puts (" x ", port);
    }
  scm_puts (">", port);
  scm_remember_upto_here_1 (state);
  return 1;                     /* non-zero means success */
}


static int vec_print (SCM vec, SCM port, scm_print_state *pstate)
{
  (void) pstate;
  const Vec c_vec = to_vec (vec);

  scm_puts ("#<vec addr: ", port);
  scm_display (from_pointer (c_vec), port);
  scm_puts (", length: ", port);
  scm_display (vec_length (vec), port);
  scm_puts (">", port);
  scm_remember_upto_here_1 (vec);
  return 1;                     /* non-zero means success */
}


static void state_init_type (void)
{
  /*
    This is the size of the struct, the actual memory pressure is much
    higher and  invisible to Guile garbage  collector. Temporary Vecs,
    matrices and such:
  */
  state_tag = scm_make_smob_type ("state", sizeof (State));
  scm_set_smob_mark (state_tag, noop_mark);
  scm_set_smob_free (state_tag, state_free);
  scm_set_smob_print (state_tag, state_print);

  /* Destroy state explicitly, when producing much garbage: */
  scm_c_define_gsubr ("state-make", 1, 0, 0, state_make);
  scm_c_define_gsubr ("state-destroy", 1, 0, 0, state_destroy);
}


static void vec_init_type (void)
{
  /* The actual memory pressure is much higher. */
  vec_tag = scm_make_smob_type ("vec", sizeof (Vec));
  scm_set_smob_mark (vec_tag, noop_mark);
  scm_set_smob_free (vec_tag, vec_free);
  scm_set_smob_print (vec_tag, vec_print);

  /* Destroy state explicitly, when producing much garbage: */
  scm_c_define_gsubr ("vec-make", 1, 0, 0, vec_make);
  scm_c_define_gsubr ("vec-make-complex", 1, 0, 0, vec_make_complex);
  scm_c_define_gsubr ("vec-destroy", 1, 0, 0, vec_destroy);

  scm_c_define_gsubr ("vec-save", 2, 0, 0, vec_save);
  scm_c_define_gsubr ("vec-load", 1, 0, 0, vec_load);
  scm_c_define_gsubr ("vec-length", 1, 0, 0, vec_length);
  scm_c_define_gsubr ("vec-ref", 2, 0, 0, vec_ref);
  scm_c_define_gsubr ("vec-set-random", 1, 0, 0, vec_set_random);
  scm_c_define_gsubr ("vec-dot", 2, 0, 0, vec_dot);
  scm_c_define_gsubr ("vec-fft", 2, 0, 0, vec_fft);
  scm_c_define_gsubr ("vec-ifft", 2, 0, 0, vec_ifft);
  scm_c_define_gsubr ("vec-fft-interp", 3, 0, 0, vec_fft_interp);
  scm_c_define_gsubr ("vec-map1", 2, 0, 0, vec_map1);
  scm_c_define_gsubr ("vec-map2", 3, 0, 0, vec_map2);
}


static SCM guile_bgy3d_solvent (SCM solvent, SCM alist)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  const ProblemData PD = problem_data (alist);

  int m;                        /* number of solvent sites */
  Site *solvent_sites;          /* solvent_sites[m] */
  char *solvent_name;

  /* Get  the  number  of   sites  and  their  parameters.   Allocates
     sol*_sites, sol*_name: */
  to_sites (solvent, &m, &solvent_sites, &solvent_name);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", solvent_name);

  /* This writes to the disk: */
  bgy3d_solve_solvent (&PD, m, solvent_sites);

  free (solvent_name);
  free (solvent_sites);

  return alist;
}


static SCM guile_hnc3d_solvent (SCM solvent, SCM settings)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  const ProblemData PD = problem_data (settings);

  int m;                        /* number of solvent sites */
  Site *solvent_sites;          /* solvent_sites[m] */
  char *solvent_name;

  /* Get  the  number  of   sites  and  their  parameters.   Allocates
     sol*_sites, sol*_name: */
  to_sites (solvent, &m, &solvent_sites, &solvent_name);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", solvent_name);

  Vec g[m][m];

  /* This writes to the disk: */
  hnc3d_solvent_solve (&PD, m, solvent_sites, g);

  free (solvent_name);
  free (solvent_sites);

  bgy3d_vec_destroy2 (m, g);

  return settings;
}


static SCM guile_bgy3d_solute (SCM solute, SCM solvent, SCM settings,
                               SCM restart) /* optional */
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  const ProblemData PD = problem_data (settings);

  int m;                        /* number of solvent sites */
  Site *solvent_sites;          /* solvent_sites[m] */
  char *solvent_name;

  int n;                        /* number of solute sites */
  Site *solute_sites;           /* solute_sites[n] */
  char *solute_name;

  /* Get  the  number  of   sites  and  their  parameters.   Allocates
     sol*_sites, sol*_name: */
  to_sites (solvent, &m, &solvent_sites, &solvent_name);
  to_sites (solute, &n, &solute_sites, &solute_name);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", solvent_name);
  PetscPrintf (PETSC_COMM_WORLD, "Solute is %s.\n", solute_name);

  /* This declares and  sets a function pointer. If  the settings dont
     specify it, it should remain NULL: */
  void (*qm_density) (int n, const real x[n][3], real rho[n]) = NULL;

  /*
    Cast is to silence the warning  here.  Note that we pass a pointer
    to  a  funptr,  void  (**)(),  as  the  function  is  supposed  to
    (eventually) set that funptr to something meaningful:
  */
  alist_getopt_funptr (settings, "qm-density", (void (**)()) &qm_density);

  /*
    A call  to bgy3d_solute_solve() takes  part of the input  from the
    disk and returns  solvent distribution in Vec g[]  (dont forget to
    destroy them).  If no additional charge distribution is associated
    with the solute  pass NULL as the function  pointer. Similarly, if
    you do not want an iterator over the solvent potential pass NULL:
  */
  Vec g[m];
  Context *medium_;

  /*
    If the argument "restart" is  not present, the caller is not going
    to  resume a  calculaiton  with with  a,  say, slightly  different
    solute  parameters.  When  the argument  is present  but  NULL the
    restart will be (eventually) used, but no restart information from
    the previous run  is available yet.  In this  case the solver puts
    such info into this  position on exit, thus, effectively returning
    one more value.
  */
  SCM next_restart;
  if (SCM_UNBNDP (restart))
    {
      bgy3d_solute_solve (&PD, m, solvent_sites, n, solute_sites, qm_density,
                          g,          /* out */
                          &medium_,   /* out, if not NULL */
                          NULL);      /* not restartable */
      next_restart = restart;         /* unspecified */
    }
  else
    {
      /*
        This is a pointer to  some structure holding restart info (ok,
        so far it is just a long Vec in disguise). This is NULL in the
        first call of a series:
      */
      Restart *restart_ = to_pointer (restart);
      bgy3d_solute_solve (&PD, m, solvent_sites, n, solute_sites, qm_density,
                          g,          /* out */
                          &medium_,   /* out */
                          &restart_); /* inout */
      next_restart = from_pointer (restart_);
    }

  free (solute_name);
  free (solute_sites);
  free (solvent_name);
  free (solvent_sites);

  /* Build a list starting from the tail: */
  SCM gs = from_vec1 (m, g);
  SCM medium = from_pointer (medium_);

  /*
    Return  multiple  values. Caller,  dont  forget  to destroy  them!
    FIXME: this is a terrible function --- the number of return values
    depends on the input!
  */
  if (SCM_UNBNDP (next_restart))
    return scm_values (scm_list_2 (gs, medium));
  else
    return scm_values (scm_list_3 (gs, medium, next_restart));
}


static SCM guile_restart_destroy (SCM restart)
{
  bgy3d_restart_destroy (to_pointer (restart));
  return from_pointer (NULL);
}


static SCM guile_hnc3d_solute (SCM solute, SCM solvent, SCM settings)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  const ProblemData PD = problem_data (settings);

  int m;                        /* number of solvent sites */
  Site *solvent_sites;          /* solvent_sites[m] */
  char *solvent_name;

  int n;                        /* number of solute sites */
  Site *solute_sites;           /* solute_sites[n] */
  char *solute_name;

  /* Get  the  number  of   sites  and  their  parameters.   Allocates
     sol*_sites, sol*_name: */
  to_sites (solvent, &m, &solvent_sites, &solvent_name);
  to_sites (solute, &n, &solute_sites, &solute_name);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", solvent_name);
  PetscPrintf (PETSC_COMM_WORLD, "Solute is %s.\n", solute_name);

  /*
    This  takes part  of  the  input from  the  disk, returns  solvent
    distribution in Vec g[] (dont forget to destroy them).:
  */
  Vec g[m];
  hnc3d_solute_solve (&PD, m, solvent_sites, n, solute_sites, g);

  free (solute_name);
  free (solute_sites);
  free (solvent_name);
  free (solvent_sites);

  /* Build a list starting from the tail: */
  SCM gs = from_vec1 (m, g);
  SCM v = from_pointer (NULL);

  /* Return multiple values. Caller, dont forget to destroy them! */
  return scm_values (scm_list_2 (gs, v));
}



static SCM guile_pot_interp (SCM iter, SCM x)
{
  double x_[1][3], v_[1];

  /* list -> array: */
  to_array (x, 3, x_[0]);

  /* printf ("guile_pot_interp: x = (% f, % f, % f)\n", */
  /*         x_[0][0], x_[0][1], x_[0][2]); */

  bgy3d_pot_interp (to_pointer (iter), 1, x_, v_);

  return scm_from_double (v_[0]);
}

static SCM guile_pot_destroy (SCM iter)
{
  bgy3d_pot_destroy (to_pointer (iter));

  return scm_from_int (0);
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
  scm_c_define_gsubr ("hnc3d-run-solvent", 2, 0, 0, guile_hnc3d_solvent);
  scm_c_define_gsubr ("hnc3d-run-solute", 3, 0, 0, guile_hnc3d_solute);
  scm_c_define_gsubr ("bgy3d-run-solvent", 2, 0, 0, guile_bgy3d_solvent);
  scm_c_define_gsubr ("bgy3d-run-solute", 3, 1, 0, guile_bgy3d_solute);
  scm_c_define_gsubr ("bgy3d-pot-interp", 2, 0, 0, guile_pot_interp);
  scm_c_define_gsubr ("bgy3d-pot-destroy", 1, 0, 0, guile_pot_destroy);
  scm_c_define_gsubr ("bgy3d-restart-destroy", 1, 0, 0, guile_restart_destroy);
  scm_c_define_gsubr ("bgy3d-rank", 0, 0, 0, guile_rank);
  scm_c_define_gsubr ("bgy3d-size", 0, 0, 0, guile_size);
  scm_c_define_gsubr ("bgy3d-test", 3, 0, 0, guile_test);

  /* Define SMOBs: */
  state_init_type();
  vec_init_type();

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
    Set global verbosity early enough.  This is the only short option!
    Use long-options prefixed by "--" as the usual convention.
  */
  verbosity = bgy3d_getopt_test ("-v"); /* extern */

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
