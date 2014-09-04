/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

#include <libguile.h>
#include <minpack.h>            /* lmdif1_() */
#include "libbgy3d.h"           /* public API */
#include "bgy3d.h"              /* private */
#include "bgy3d-getopt.h"       /* Implementation here */
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* bgy3d_solvent_get() */
#include "bgy3d-pure.h"
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-impure.h"       /* bgy3d_solve_with_solute */
#include "hnc3d.h"              /* hnc3d_solute_solve() */
#include "bgy3d-vec.h"          /* bgy3d_vec_save, bgy3d_vec_load */
#include "bgy3d-fft.h"          /* bgy3d_fft_test() */
#include "bgy3d-fftw.h"         /* bgy3d_fft_interp() */
#include "rism-dst.h"           /* rism_dst() */
#include "rism-rdf.h"           /* rism_rdf() */
#include "rism.h"               /* rism_solvent() */
#include "eos.h"                /* eos_alj(), etc. */
#include "lebed/lebed.h"        /* genpts() */
#include "bgy3d-guile.h"

#ifdef WITH_FFTW_THREADS
#include <fftw3.h>              /* fftw_init_threads(), fftw_cleanup_threads */

/* If non-zero,  tell FFTW  to use that  many threads. Use  stock FFTW
   otherwise: */
static int nthreads = WITH_FFTW_THREADS;
#endif



static char helptext[] = "BGY3d Guile.\n";

/*
  C-pointers are handled  on the Scheme side as  "pointer objects" not
  just  plain integers.  These  two convert  between two  worlds.  The
  const qualifier might  be misleading, as we cannot  be sure what the
  caller does with the result:
*/
static SCM
from_pointer (const void *ptr)
{
  return scm_from_pointer ((void *) ptr, NULL);
}


/* You cannot pass 0 here, use %null-pointer instead: */
static void*
to_pointer (SCM ptr)
{
  /* FIXME: is scm_to_pointer() available yet? */
  if (!SCM_POINTER_P (ptr))
    scm_misc_error              /* longjmp! */
      (__func__, "not a pointer object ~A", scm_list_1 (ptr));

  return SCM_POINTER_VALUE (ptr);
}


/* True if the key is present: */
static bool
alist_getopt_test (SCM alist, const char *key)
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  return scm_is_true (kv);
}

/* Update SCM value  with the entry from an  association list or leave
   it  unchanged if  there is  no corresponding  entry.  Alist  is not
   modified. */
static bool
alist_getopt_scm (SCM alist, const char *key, SCM *val)
{
  SCM kv = scm_assoc (scm_from_locale_symbol (key), alist);

  bool test = scm_is_true (kv);

  if (test)
    *val = scm_cdr (kv);

  return test;
}

/* Update a bool  with the entry from an association  list or leave it
   unchanged  if  there  is  no  corresponding entry.   Alist  is  not
   modified. */
static bool
alist_getopt_bool (SCM alist, const char *key, bool *val)
{
  SCM scm;
  bool test = alist_getopt_scm (alist, key, &scm);

  if (test)
    *val = scm_to_bool (scm);

  return test;
}


/* Update an integer with the  entry from an association list or leave
   it  unchanged if  there is  no  corresponding entry.  Alist is  not
   modified. */
static bool
alist_getopt_int (SCM alist, const char *key, int *val)
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
static bool
alist_getopt_real (SCM alist, const char *key, double *val)
{
  SCM scm;
  bool test = alist_getopt_scm (alist, key, &scm);

  if (test)
    *val = scm_to_double (scm);

  return test;
}

/* Update a function  pointer with the entry from  an association list
   or leave it unchanged if there is no corresponding entry: */
static bool
alist_getopt_funptr (SCM alist, /* intent(in) */
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

/* Lookup a name in a module: */
static SCM
lookup (const char *module, const char *name)
{
  SCM mod = scm_c_resolve_module (module);
  return scm_variable_ref (scm_c_module_lookup (mod, name));
}


/* Lookup dynvar: */
static SCM
guile_get_settings ()
{
  return scm_fluid_ref (lookup ("guile bgy3d", "*settings*"));
}


bool
bgy3d_getopt_test (const char key[])
{
  return alist_getopt_test (guile_get_settings (), key);
}


bool
bgy3d_getopt_bool (const char key[], bool *val)
{
  return alist_getopt_bool (guile_get_settings (), key, val);
}


bool
bgy3d_getopt_scm (const char key[], SCM *val)
{
  return alist_getopt_scm (guile_get_settings (), key, val);
}


bool
bgy3d_getopt_int (const char key[], int *val)
{
  return alist_getopt_int (guile_get_settings (), key, val);
}


bool
bgy3d_getopt_real (const char key[], double *val)
{
  return alist_getopt_real (guile_get_settings (), key, val);
}


bool
bgy3d_getopt_string (const char key[], int len, char val[len])
{
  SCM str;
  bool ok = alist_getopt_scm (guile_get_settings (), key, &str);

  if (ok)
    {
      /* Does not null-terminate, returns the number of bytes required
         for the string not counting the terminating \0: */
      int n = scm_to_locale_stringbuf (str, val, len);

      if (n < len)
        val[n] = '\0';
      else
        val[len - 1] = '\0'; /* FIXME: truncation here! */
    }

  return ok;
}


/*
  Get  problem data  from the  association  list (set  e.g.  from  the
  command line).  This interface evolved from asking PETSC environment
  to  consulting  the settings  fluid  passed  here  explicitly as  an
  association list over  time. It is not always  convenient to restart
  executable to change a parameter.

  Do not forget to set the defaults, the environment may not contain a
  setting for every field!  FIXME:  some parameters are handled by the
  working code itself having potentially its own defaults.

  Will longjmp() on error!
*/
static ProblemData
problem_data (SCM alist)
{
  ProblemData PD;

  /* Inverse temperature: */
  PD.beta = 1.6889;

  /* Density: */
  PD.rho = 0.3;

  /* Mixing parameter: */
  PD.lambda = 0.1;

  /*
    Initial  interaction scaling  factor.  The code  used to  automate
    "blending-in"  some interaction by  iterating scaling  factor from
    some low initial value to 1.  Default was 0.0 originally, but that
    makes necessary to always specify --damp-start when 1.0 you do not
    want blending.  I think 1.0 is a more reasonable default:
  */
  PD.damp = 1.0;

  /* Number of total iterations: */
  PD.max_iter = 100;

  /* Norm tolerance for convergence test: */
  PD.norm_tol = 1.0e-12;

  /* Overwrite  defaults with  the  data provided  in the  association
     list. First real (double precision) options: */
  alist_getopt_real (alist, "rho", &PD.rho);
  alist_getopt_real (alist, "beta", &PD.beta);
  alist_getopt_real (alist, "norm-tol", &PD.norm_tol);
  alist_getopt_real (alist, "lambda", &PD.lambda);
  alist_getopt_real (alist, "damp-start", &PD.damp);

  assert (PD.lambda >= 0.0);
  assert (PD.lambda <= 1.0);

  /* (half of the) box size: */
  real length = 12.0;
  alist_getopt_real (alist, "L", &length);

  for (int i = 0; i < 3; i++)
    PD.L[i] = 2 * length;

  /* Grid points in 1 dimension */
  int n = 32;
  alist_getopt_int (alist, "N", &n);
  assert (n > 0);

  for (int i = 0; i < 3; i++)
    PD.N[i] = n;

  alist_getopt_int (alist, "max-iter", &PD.max_iter);

  /* Preserve invariants or get rid of redundancies: */
  for (int i = 0; i < 3; i++)
    PD.h[i] = PD.L[i] / PD.N[i];

  /*
    At this point N, h and L have consistent values.  FIXME: N^2 + N^2
    + N^2 should not overflow  in 3D, this condition ensures that 3N^2
    <  2^31.  Also  N^3  should not  ne  too large.   But  for 1D  the
    restriction is void so do not abort if N^3 >= 2^31:
  */
  assert (n < 26755);
  // assert (n < 1291);

  /*
    Radial  parameters for use  in 1D  RISM, were  previousely derived
    like this  directly from L[] and  N[]. Note that the  3D code when
    calling 1D RISM will scale rmax and nrad up, by default.
  */
  PD.rmax = MAX (MAX (PD.L[0], PD.L[1]), PD.L[2]) / 2;
  PD.nrad = MAX (MAX (PD.N[0], PD.N[1]), PD.N[2]);
  alist_getopt_real (alist, "rmax", &PD.rmax);
  alist_getopt_int (alist, "nrad", &PD.nrad);

  /*
    Closure is only  used in HNC-like methods. Supply  default for the
    case the settings do not provide it.

    Closure is  a symbol,  in capital letters.   Only do  something if
    closure is specified in the association list, otherwise assume the
    default:
  */
  PD.closure = CLOSURE_HNC;
  SCM closure;
  if (alist_getopt_scm (alist, "closure", &closure))
    {
      const SCM hnc = scm_from_locale_symbol ("HNC");
      const SCM kh = scm_from_locale_symbol ("KH");
      const SCM py = scm_from_locale_symbol ("PY");
      const SCM pse0 = scm_from_locale_symbol ("PSE0");
      const SCM pse1 = scm_from_locale_symbol ("PSE1");
      const SCM pse2 = scm_from_locale_symbol ("PSE2");
      const SCM pse3 = scm_from_locale_symbol ("PSE3");
      const SCM pse4 = scm_from_locale_symbol ("PSE4");
      const SCM pse5 = scm_from_locale_symbol ("PSE5");
      const SCM pse6 = scm_from_locale_symbol ("PSE6");
      const SCM pse7 = scm_from_locale_symbol ("PSE7");

      if (scm_is_true (scm_equal_p (closure, hnc)))
        PD.closure = CLOSURE_HNC;
      else if (scm_is_true (scm_equal_p (closure, kh)))
        PD.closure = CLOSURE_KH;
      else if (scm_is_true (scm_equal_p (closure, py)))
        PD.closure = CLOSURE_PY;
      else if (scm_is_true (scm_equal_p (closure, pse0)))
        PD.closure = CLOSURE_PSE0;
      else if (scm_is_true (scm_equal_p (closure, pse1)))
        PD.closure = CLOSURE_PSE1;
      else if (scm_is_true (scm_equal_p (closure, pse2)))
        PD.closure = CLOSURE_PSE2;
      else if (scm_is_true (scm_equal_p (closure, pse3)))
        PD.closure = CLOSURE_PSE3;
      else if (scm_is_true (scm_equal_p (closure, pse4)))
        PD.closure = CLOSURE_PSE4;
      else if (scm_is_true (scm_equal_p (closure, pse5)))
        PD.closure = CLOSURE_PSE5;
      else if (scm_is_true (scm_equal_p (closure, pse6)))
        PD.closure = CLOSURE_PSE6;
      else if (scm_is_true (scm_equal_p (closure, pse7)))
        PD.closure = CLOSURE_PSE7;
      else
        scm_misc_error          /* longjmp! */
          (__func__, "no such closure ~A", scm_list_1 (closure));
    }

  return PD;
}

/* FIXME: declaration in bgy3d.h: */
ProblemData
bgy3d_problem_data (void)
{
  return problem_data (guile_get_settings ());
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
  to_double1 (x, 3, S.x);

  S.sigma = scm_to_double (scm_car (ff));
  S.epsilon = scm_to_double (scm_cadr (ff));
  S.charge = scm_to_double (scm_caddr (ff));
  S.site = s;                  /* Link to Scheme representation too */

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
static scm_t_bits guile_state_tag;
static scm_t_bits guile_vec_tag;


static SCM from_state (const State *BHD)
{
  SCM obj;
  SCM_NEWSMOB (obj, guile_state_tag, BHD);
  return obj;
}


static SCM from_vec (Vec vec)
{
  SCM obj;
  SCM_NEWSMOB (obj, guile_vec_tag, vec);
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
  scm_assert_smob_type (guile_state_tag, state);

  return (State*) SCM_SMOB_DATA (state);
}


static Vec to_vec (SCM vec)
{
  scm_assert_smob_type (guile_vec_tag, vec);

  return (Vec) SCM_SMOB_DATA (vec);
}


static SCM guile_state_make (SCM alist)
{
  /* This sets defaults, eventually modified from the command line and
     updated by the entries from the association list: */
  ProblemData *PD = malloc (sizeof *PD); /* free() in guile_state_free() */
  *PD = problem_data (alist);            /* FIXME: leak on longjmp! */

  return from_state (bgy3d_state_make (PD));
}


static SCM guile_vec_make (SCM state)
{
  State *BHD = to_state (state);
  Vec vec = vec_create (BHD->da);
  return from_vec (vec);
}


static SCM guile_vec_make_complex (SCM state)
{
  State *BHD = to_state (state);
  Vec vec = vec_create (BHD->dc);
  return from_vec (vec);
}


static size_t guile_state_free (SCM state)
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


static size_t guile_vec_free (SCM vec)
{
  Vec c_vec = to_vec (vec);
  if (c_vec)
    vec_destroy (&c_vec);
  return 0;
}


static SCM guile_state_destroy (SCM state)
{
  assert (to_state (state) != NULL);
  guile_state_free (state);
  SCM_SET_SMOB_DATA (state, NULL);
  return SCM_UNSPECIFIED;
}


static SCM guile_vec_destroy (SCM vec)
{
  assert (to_vec (vec) != NULL);
  guile_vec_free (vec);
  SCM_SET_SMOB_DATA (vec, NULL);
  return SCM_UNSPECIFIED;
}


static SCM noop_mark (SCM smob)
{
  /* We dont store SCM values in smobs yet: */
  (void) smob;
  return SCM_BOOL_F;
}


static SCM guile_vec_save (SCM path, SCM vec)
{
  char *c_path = scm_to_locale_string (path); /* free() it! */

  bgy3d_vec_save (c_path, to_vec (vec));

  free (c_path);

  return SCM_UNDEFINED;
}


static SCM guile_vec_load (SCM path)
{
  char *c_path = scm_to_locale_string (path); /* free() it! */

  Vec vec = bgy3d_vec_load (c_path);

  free (c_path);

  return from_vec (vec);
}


static SCM guile_vec_length (SCM vec)
{
  return scm_from_int (vec_size (to_vec (vec)));
}

/* An inefficient way of getting just one value, vec[ix], even if that
   value is not stored locally. Collective. */
static SCM guile_vec_ref (SCM vec, SCM ix)
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
  comm_allreduce (1, vals);

  return scm_from_double (vals[0]);
}


/* Desctructively updates Vec */
static SCM guile_vec_shift_x (SCM vec, SCM shift)
{
  VecShift (to_vec (vec), scm_to_double (shift));

  return vec;
}


/* Desctructively updates Vec */
static SCM guile_vec_scale_x (SCM vec, SCM scale)
{
  VecScale (to_vec (vec), scm_to_double (scale));

  return vec;
}


static SCM guile_vec_dot (SCM x, SCM y)
{
  return scm_from_double (vec_dot (to_vec (x), to_vec (y)));
}


static SCM guile_vec_fft (SCM state, SCM x)
{
  State *BHD = to_state (state);
  SCM y = guile_vec_make_complex (state);
  Vec y_ = to_vec (y);

  const int N3 = BHD->PD->N[0] * BHD->PD->N[1] * BHD->PD->N[2];

  MatMult (BHD->fft_mat, to_vec (x), y_);
  VecScale (y_, 1.0 / sqrt (N3));

  return y;
}


static SCM guile_vec_ifft (SCM state, SCM y)
{
  State *BHD = to_state (state);
  SCM x = guile_vec_make (state);
  Vec x_ = to_vec (x);

  const int N3 = BHD->PD->N[0] * BHD->PD->N[1] * BHD->PD->N[2];

  MatMultTranspose (BHD->fft_mat, to_vec (y), x_);
  VecScale (x_, 1.0 / sqrt (N3));

  return x;
}


static SCM guile_vec_fft_interp (SCM state, SCM Y, SCM x)
{
  State *BHD = to_state (state);
  Vec Y_ = to_vec (Y);
  const int N3 = BHD->PD->N[0] * BHD->PD->N[1] * BHD->PD->N[2];

  double x_[1][3], y[1];
  for (int i = 0; i < 3; i ++)
    {
      x_[0][i] = scm_to_double (scm_car (x));
      x = scm_cdr (x);
    }

  bgy3d_fft_interp (BHD->fft_mat, Y_, 1, x_, y);

  /* FIXME: FFT in BGY3d code is unnormalized: */
  return scm_from_double (y[0] * sqrt (N3));
}


/* (rism-rdf domain center radial-mesh angular-order) */
static SCM
guile_rism_rdf (SCM state, SCM g, SCM center, SCM radial_mesh, SCM angular_order)
{
  /* RDF center: */
  double a[3];
  to_double1 (center, 3, a);

  /* Number of radial pioints: */
  const int n = scm_to_int (scm_length (radial_mesh));

  /* Radial points: */
  double r[n];
  to_double1 (radial_mesh, n, r);

  const int m = scm_to_int (angular_order);

  /* Output array: */
  double rdf[n];

  /* Does the real work: */
  rism_rdf (to_state (state), to_vec (g), a, n, r, m, rdf);

  return from_double1 (n, rdf);
}


/* (vec-moments domain vec) */
static SCM
guile_vec_moments (SCM state, SCM g)
{
  real m0, m1[3], m2[3][3];

  /* Does the real work: */
  bgy3d_moments (to_state (state), to_vec (g), &m0, m1, m2);

  return scm_list_3 (scm_from_double (m0),
                     from_double1 (3, m1),
                     from_double2 (3, 3, m2));
}


static SCM
guile_genpts (SCM M)
{
  const int m = scm_to_int (M);
  double x[m][3], w[m];

  const int n = genpts (m, x, w);
  assert (n <= m);

  /* n == m */
  return scm_values (scm_list_2 (from_double2 (n, 3, x), from_double1 (n, w)));
}


static SCM guile_vec_set_random (SCM x)
{
  VecSetRandom (to_vec (x), NULL);

  return x;
}

/*
  NOTE:  scm_from_double()   allocates  on  heap   producing  lots  of
  garbage. A  good thing  about it is  that it  also causes GC  to run
  regularly and collect Vec smobs.
*/
static SCM guile_vec_map1 (SCM f, SCM x)
{
  Vec x_ = to_vec (x);
  Vec y_ = vec_duplicate (x_);

  void f_ (int n, double x[n], double y[n])
  {
    for (int i = 0; i < n; i++)
      y[i] = scm_to_double (scm_call_1 (f, scm_from_double (x[i])));
  }
  vec_app2 (f_, x_, y_);

  return from_vec (y_);
}


static SCM guile_vec_map2 (SCM f, SCM x, SCM y)
{
  Vec x_ = to_vec (x);
  Vec y_ = to_vec (y);
  Vec z_ = vec_duplicate (x_);

  void f_ (int n, double x[n], double y[n], double z[n])
  {
    for (int i = 0; i < n; i++)
      z[i] = scm_to_double
        (scm_call_2 (f, scm_from_double (x[i]), scm_from_double (y[i])));
  }
  vec_app3 (f_, x_, y_, z_);

  return from_vec (z_);
}


static int guile_state_print (SCM state, SCM port, scm_print_state *pstate)
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


static int guile_vec_print (SCM vec, SCM port, scm_print_state *pstate)
{
  (void) pstate;
  const Vec c_vec = to_vec (vec);

  scm_puts ("#<vec addr: ", port);
  scm_display (from_pointer (c_vec), port);
  /* Someone  will try to  display a  Vec that  has been  destroyed in
     Scheme: */
  if (c_vec)
    {
      scm_puts (", length: ", port);
      scm_display (guile_vec_length (vec), port);
    }
  scm_puts (">", port);
  scm_remember_upto_here_1 (vec);
  return 1;                     /* non-zero means success */
}

#if 1   /* FIXME: move them out of the way into a separate file ... */
static inline SCM app1 (void (*f)(size_t n, double y[n], const double x[n]),
                        SCM x)
{
  scm_t_array_handle hx, hy;
  size_t nx, ny;
  ssize_t dx, dy;

  const double *x_ = scm_uniform_vector_elements (x, &hx, &nx, &dx);
  assert (dx == 1);

  SCM y = scm_make_f64vector (scm_from_int (nx), SCM_UNDEFINED);

  double *y_ = scm_uniform_vector_writable_elements (y, &hy, &ny, &dy);
  assert (dy == 1);
  assert (ny == nx);

  f (nx, y_, x_);

  scm_array_handle_release (&hy);
  scm_array_handle_release (&hx);

  return y;
}


static inline SCM f64_map1 (double (*f)(double), SCM x)
{
  void g (size_t n, double y[n], const double x[n])
  {
    for (size_t i = 0; i < n; i++)
      y[i] = f (x[i]);
  }
  return app1 (g, x);
}


static inline SCM f64_map2 (double (*f)(double, double), SCM x, SCM y)
{
  scm_t_array_handle hx, hy, hz;
  size_t nx, ny, nz;
  ssize_t dx, dy, dz;

  const double *x_ = scm_uniform_vector_elements (x, &hx, &nx, &dx);
  assert (dx == 1);

  const double *y_ = scm_uniform_vector_elements (y, &hy, &ny, &dy);
  assert (dy == 1);
  assert (ny == nx);

  SCM z = scm_make_f64vector (scm_from_int (nx), SCM_UNDEFINED);

  double *z_ = scm_uniform_vector_writable_elements (z, &hz, &nz, &dz);
  assert (dz == 1);
  assert (nz == nx);

  for (size_t i = 0; i < nz; i++)
    z_[i] = f (x_[i], y_[i]);

  scm_array_handle_release (&hz);
  scm_array_handle_release (&hy);
  scm_array_handle_release (&hx);

  return z;
}


static SCM guile_f64_scale (SCM a, SCM x)
{
  const double a_ = scm_to_double (a);
  double f (double x)
  {
    return a_ * x;
  }
  return f64_map1 (f, x);
}


static SCM guile_f64_add (SCM x, SCM y)
{
  double f (double x, double y)
  {
    return x + y;
  }
  return f64_map2 (f, x, y);
}


static SCM guile_f64_mul (SCM x, SCM y)
{
  double f (double x, double y)
  {
    return x * y;
  }
  return f64_map2 (f, x, y);
}


static SCM guile_f64_dst (SCM x)
{
  return app1 (rism_dst, x);
}
#endif  /* ... staff that belongs elsewhere. */


#define EXPORT(name, req, opt, rst, func)               \
  (scm_c_define_gsubr (name, req, opt, rst, func),      \
   scm_c_export (name, NULL))

static void guile_init_state_type (void)
{
  /*
    This is the size of the struct, the actual memory pressure is much
    higher and  invisible to Guile garbage  collector. Temporary Vecs,
    matrices and such:
  */
  guile_state_tag = scm_make_smob_type ("state", sizeof (State));
  scm_set_smob_mark (guile_state_tag, noop_mark);
  scm_set_smob_free (guile_state_tag, guile_state_free);
  scm_set_smob_print (guile_state_tag, guile_state_print);

  /* Destroy state explicitly, when producing much garbage: */
  EXPORT ("state-make", 1, 0, 0, guile_state_make);
  EXPORT ("state-destroy", 1, 0, 0, guile_state_destroy);
}


static void guile_init_vec_type (void)
{
  /* The actual memory pressure is much higher. */
  guile_vec_tag = scm_make_smob_type ("vec", sizeof (Vec));
  scm_set_smob_mark (guile_vec_tag, noop_mark);
  scm_set_smob_free (guile_vec_tag, guile_vec_free);
  scm_set_smob_print (guile_vec_tag, guile_vec_print);

  /* Destroy state explicitly, when producing much garbage: */
  EXPORT ("vec-make", 1, 0, 0, guile_vec_make);
  EXPORT ("vec-make-complex", 1, 0, 0, guile_vec_make_complex);
  EXPORT ("vec-destroy", 1, 0, 0, guile_vec_destroy);

  EXPORT ("vec-save", 2, 0, 0, guile_vec_save);
  EXPORT ("vec-load", 1, 0, 0, guile_vec_load);
  EXPORT ("vec-length", 1, 0, 0, guile_vec_length);
  EXPORT ("vec-ref", 2, 0, 0, guile_vec_ref);
  EXPORT ("vec-set-random", 1, 0, 0, guile_vec_set_random);
  EXPORT ("vec-dot", 2, 0, 0, guile_vec_dot);
  EXPORT ("vec-fft", 2, 0, 0, guile_vec_fft);
  EXPORT ("vec-ifft", 2, 0, 0, guile_vec_ifft);
  EXPORT ("vec-fft-interp", 3, 0, 0, guile_vec_fft_interp);
  EXPORT ("vec-map1", 2, 0, 0, guile_vec_map1);
  EXPORT ("vec-map2", 3, 0, 0, guile_vec_map2);
  EXPORT ("vec-moments", 2, 0, 0, guile_vec_moments);
  EXPORT ("vec-shift!", 2, 0, 0, guile_vec_shift_x);
  EXPORT ("vec-scale!", 2, 0, 0, guile_vec_scale_x);

  EXPORT ("rism-rdf", 5, 0, 0, guile_rism_rdf);
  EXPORT ("genpts", 1, 0, 0, guile_genpts);
  EXPORT ("f64dst", 1, 0, 0, guile_f64_dst);
  EXPORT ("f64+", 2, 0, 0, guile_f64_add);
  EXPORT ("f64*", 2, 0, 0, guile_f64_mul);
  EXPORT ("f64scale", 2, 0, 0, guile_f64_scale);
  EXPORT ("eos-alj", 2, 0, 0, eos_alj);
  EXPORT ("eos-alj-res", 2, 0, 0, eos_alj_res);
  EXPORT ("eos-plj", 2, 0, 0, eos_plj);
  EXPORT ("eos-ulj", 2, 0, 0, eos_ulj);
}


/* E.g. hnc3d_solvent_solve() or bgy3d_solvent_solve(): */
typedef
void (*SV) (const ProblemData *PD, int m, const Site solvent[m], Vec g[m][m]);

/* Decode SCM input,  encode output to SCM. The  first argument is the
   actual solver that operates with C-types. */
static SCM
run_solvent (SV solve_solvent, SCM solvent)
{
  /* Lookup dynvar: */
  SCM settings = guile_get_settings ();

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
  PRINTF ("Solvent is %s.\n", solvent_name);

  local Vec g[m][m];

  /* This writes to the disk: */
  solve_solvent (&PD, m, solvent_sites, g);

  free (solvent_name);
  free (solvent_sites);

  vec_destroy2 (m, g);

  return settings;
}


static SCM
guile_bgy3d_solvent (SCM solvent)
{
  return run_solvent (bgy3d_solve_solvent, solvent);
}


static SCM
guile_hnc3d_solvent (SCM solvent)
{
  return run_solvent (hnc3d_solvent_solve, solvent);
}


/* E.g. hnc3d_solute_solve() or bgy3d_solute_solve(): */
typedef
void (*SU) (const ProblemData *PD,
            const int m, const Site solvent[m],
            const int n, const Site solute[n],
            void (*density)(int k, const real x[k][3], real rho[k]),
            SCM *dict,          /* inout, alist */
            Vec g[m],           /* out */
            Context **medium,   /* out */
            Restart **restart); /* inout */



/* Decode SCM input,  encode output to SCM. The  first argument is the
   actual solver that operates with C-types. */
static SCM
run_solute (SU solute_solve, SCM solute, SCM solvent, SCM restart)
{
  /* Lookup dynvar: */
  SCM settings = guile_get_settings ();

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
  PRINTF ("Solvent is %s.\n", solvent_name);
  PRINTF ("Solute is %s.\n", solute_name);

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
    A call to solute_solve() takes part of the input from the disk and
    returns solvent  distribution in Vec  g[] (dont forget  to destroy
    them).  If  no additional  charge distribution is  associated with
    the solute pass NULL as the function pointer. Similarly, if you do
    not want an iterator over the solvent potential pass NULL:
  */
  Vec g[m];
  Context *medium_;

  /*
    If the  argument SCM  restart is SCM  NULL no  restart information
    from the previous run is available yet.  This is a pointer to some
    structure holding restart  info (ok, so far it is  just a long Vec
    in disguise). This is NULL in the first call of a series:
  */
  Restart *restart_ = to_pointer (restart); /* maybe NULL */

  /* The code will fill the dictionary with results: */
  SCM dict = SCM_EOL;

  solute_solve (&PD, m, solvent_sites, n, solute_sites, qm_density,
                &dict,          /* inout, alist */
                g,              /* out */
                &medium_,       /* out */
                &restart_);     /* inout */

  /* Not NULL if solver supports restarting: */
  restart = from_pointer (restart_);

  free (solute_name);
  free (solute_sites);
  free (solvent_name);
  free (solvent_sites);

  /* Build a list starting from the tail: */
  SCM gs = from_vec1 (m, g);
  SCM medium = from_pointer (medium_);

  /* Extend  result   dictionary.   Caller  is   supposed  to  destroy
     these! */
  dict = scm_acons (scm_from_locale_symbol ("GUV"), gs, dict);
  dict = scm_acons (scm_from_locale_symbol ("POTENTIAL"), medium, dict);
  dict = scm_acons (scm_from_locale_symbol ("RESTART"), restart, dict);

  return dict;
}


static SCM
guile_bgy3d_solute (SCM solute, SCM solvent, SCM restart)
{
  return run_solute (bgy3d_solute_solve, solute, solvent, restart);
}


static SCM
guile_hnc3d_solute (SCM solute, SCM solvent, SCM restart)
{
  return run_solute (hnc3d_solute_solve, solute, solvent, restart);
}


static SCM guile_restart_destroy (SCM restart)
{
  bgy3d_restart_destroy (to_pointer (restart));
  return from_pointer (NULL);
}


static SCM guile_rism_solvent (SCM solvent)
{
  /* Lookup dynvar: */
  SCM settings = guile_get_settings ();

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
  if (verbosity > 0)
    PRINTF (" # Solvent is %s.\n", solvent_name);


  /* Always use this function to derive number of radial points: */
  const int nrad = PD.nrad;

  /* Guile m x m x nrad array of doubles: */
  SCM chi_fft = scm_make_typed_array (scm_from_locale_symbol ("f64"),
                                      SCM_UNSPECIFIED,
                                      scm_list_3 (scm_from_int (m),
                                                  scm_from_int (m),
                                                  scm_from_int (nrad)));

  /* This association list will contain essential results: */
  SCM retval;
  {
    scm_t_array_handle handle;
    scm_array_get_handle (chi_fft, &handle);

    /*
      This  should have  as much  space as  a real  array  declared as
      double x_buf[m][m][nrad].  Void* is to silence the warnings:
    */
    void *x_buf = scm_array_handle_f64_writable_elements (&handle);

    /*
      Actual  pure  solvent   calculation  here.   NULL  indicates  an
      optional output  argument. Supply NULL  if you dont  need either
      solvent indirect correlation or solvent susceptibility:
    */
    rism_solvent (&PD, m, solvent_sites, NULL, x_buf, &retval);

    scm_array_handle_release (&handle);
  }

  free (solvent_name);
  free (solvent_sites);

  /* Cons a  new entry onto  the association list. Note  that printing
     this m x m x nrad array may take a lot of screen estate: */
  retval = scm_acons (scm_from_locale_symbol ("susceptibility"),
                      chi_fft,
                      retval);

  return retval;
}



static SCM guile_rism_solute (SCM solute, SCM solvent, SCM chi_fft)
{
  /* Lookup dynvar: */
  SCM settings = guile_get_settings ();

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
  if (verbosity > 0)
    {
      PRINTF (" # Solvent is %s.\n", solvent_name);
      PRINTF (" # Solute is %s.\n", solute_name);
    }

  /* This association list will contain essential results: */
  SCM retval;
  {
    scm_t_array_handle handle;

    /*
      By default, there is no solvent susceptibility from a prior pure
      solvent  calculation  to  be  used here.  The  Fortran  function
      rism_solute() will have to run a pure solvent calculation too in
      this case:
    */
    const real *x_buf = NULL;

    /*
      Only if the caller supplied a suitable Guile array, then we pass
      its contents  further.  The caller  is responsible to  make sure
      that susceptibility corresponds to  the actual solvent and other
      settings (notably radial dimensions):
    */
    if (!SCM_UNBNDP (chi_fft))
      {
        scm_array_get_handle (chi_fft, &handle);

        /*
          This should have  as much space as a  real array declared as
          double x_buf[m][m][nrad].  The buffer is  read-only here, so
          we  cast the  pointer  to (void*)  when  passing further  to
          silence the warning.
        */
        x_buf = scm_array_handle_f64_elements (&handle);
      }

    /* Actual solute/solvent calculation here: */
    rism_solute (&PD, n, solute_sites, m, solvent_sites,
                 (void*) x_buf, &retval);

    if (!SCM_UNBNDP (chi_fft))
      scm_array_handle_release (&handle);
  }


  free (solvent_name);
  free (solvent_sites);

  free (solute_name);
  free (solute_sites);

  return retval;
}



static SCM
guile_rism_self_energy (SCM solute, SCM species)
{
  int n;                        /* number of solute sites */
  Site *solute_sites;           /* solute_sites[n] */
  char *solute_name;

  to_sites (solute, &n, &solute_sites, &solute_name);

  int spec[n];

  /* SCM species is a list of integers: */
  to_int1 (species, n, spec);

  /* Multiple values: */
  SCM eg = rism_self_energy (n, solute_sites, spec);

  free (solute_name);
  free (solute_sites);

  return eg;
}



static SCM guile_pot_interp (SCM iter, SCM x)
{
  double x_[1][3], v_[1];

  /* list -> array: */
  to_double1 (x, 3, x_[0]);

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

/* Return MPI runk in comm_world_petsc: */
static SCM guile_comm_rank (void)
{
  return scm_from_int (comm_rank ());
}


/* Return MPI size of comm_world_petsc: */
static SCM guile_comm_size (void)
{
  return scm_from_int (comm_size ());
}

static SCM
guile_comm_set_parallel_x (SCM flag)
{
  return scm_from_bool (comm_set_parallel_x (scm_to_bool (flag)));
}


#if SCM_MAJOR_VERSION > 1       /* Bytevectors not availbale in 1.8 */
/* read!  and write!  for use in make-custom-binary-input- and
   output-port. */
static SCM
guile_comm_bcast_x (SCM root, SCM bv, SCM start, SCM len)
{
  const size_t i = scm_to_size_t (start);
  const size_t n = scm_to_size_t (len);

  assert (scm_is_bytevector (bv));
  const size_t size = SCM_BYTEVECTOR_LENGTH (bv);
  char *buf = (void*) SCM_BYTEVECTOR_CONTENTS (bv);
  assert (i + n <= size);

  /* Extract MPI_Comm, verifies the type: */
  const MPI_Comm comm = comm_world_petsc; // scm_to_comm (world);

  int rank, ierr;
  ierr = MPI_Comm_rank (comm, &rank);
  assert (MPI_SUCCESS==ierr);

  // printf ("ZZZ: rank=%d, n=%zu, size=%zu\n", rank, n, size);

  const int root_rank = scm_to_int (root);

  /* size_t varies between 32/64 platforms, set MPI_SIZE_T here: */
  MPI_Datatype MPI_SIZE_T;
  if (sizeof (size_t) == sizeof (unsigned long))
    MPI_SIZE_T = MPI_UNSIGNED_LONG;
  else if (sizeof (size_t) == sizeof (unsigned int))
    MPI_SIZE_T = MPI_UNSIGNED;
  else
    assert (0);

  /*
    Agree  the  size of  the  data that  can  be  transferred in  this
    transaction.   There  are no  guarantees  about  the  size of  the
    read/write buffers, so taking the minimum is the best we can do:
  */
  size_t nn = n;
  ierr = MPI_Allreduce (MPI_IN_PLACE, &nn, 1, MPI_SIZE_T, MPI_MIN, comm);
  assert (MPI_SUCCESS==ierr);

  /*
    We could have done MPI_Bcast (&nn, 1, MPI_SIZE_T, root_rank, comm)
    instead, but we will have a  problem if the reader wants less than
    the writer  side has to offer.   In this case the  reader does not
    supply enough  space to  receive broadcasted data.   The assertion
    below fires indeed. That is why we need to allreduce/min the input
    n over  all workers. Then  the assertion holds trivially.   On the
    other hand  the writer is  not guaranteed to  be able to  sink all
    data in one transaction.

    The write size is not  constrained from above as a (put-bytevector
    port BV) will pass BV  through, irrespective of its size. It seems
    there is  no buffering  on the  writer side as  of 2.0.5.  Also, a
    simple  (write data  port) will  end calling  this  procedure many
    times with just a few bytes.
  */
  assert (nn <= n);

  ierr = MPI_Bcast (&buf[i], nn, MPI_CHAR, root_rank, comm);
  assert (MPI_SUCCESS==ierr);

  /* If nn  < n on the  writer side, the caller  will initiate another
     transaction: */
  return scm_from_size_t (nn);
}
#endif


/*
  Find a list of  numbers x such that the sum of  squares of f(x) list
  is minmal. Start with x = x0.

  This is the signature of MIPACK procedure:

    SUBROUTINE LMDIF1 (FCN, M, N, X, FVEC, TOL, INFO, IWA, WA, LWA)
    INTEGER M, N, INFO, LWA
    INTEGER IWA(N)
    DOUBLE PRECISION TOL
    DOUBLE PRECISION X(N), FVEC(M), WA(LWA)

  where

    SUBROUTINE FCN (M, N, X, FVEC, IFLAG)
    INTEGER M, N, IFLAG
    DOUBLE PRECISION X(N), FVEC(M)
*/
static SCM
guile_least_squares (SCM f, SCM x0)
{
  int n = scm_to_int (scm_length (x0));
  SCM f0 = scm_call_1 (f, x0);  /* FIXME: need to know the shape */
  int m = scm_to_int (scm_length (f0));

  void fcn (int *m, int *n, double x[], double fvec[], int *iflag)
  {
    (void) iflag;
    SCM x1 = from_double1 (*n, x);
    SCM f1 = scm_call_1 (f, x1);
    to_double1 (f1, *m, fvec);
  }

  double x[n];
  to_double1 (x0, n, x);

  double fvec[m];
  double tol = 1.0e-7;          /* FIXME: literal here! */
  int info, iwa[n];
  int lwa = m * n + 5 * n + m;
  double wa[lwa];

  /* The joys of calling F77 ... */
  lmdif1_ (fcn, &m, &n, x, fvec, &tol, &info, iwa, wa, &lwa);

  assert (info == 1 || info == 2 || info == 3); /* FXIME: 4 too?*/
  /*
    INFO = 0 Improper input parameters.

    INFO = 1 Algorithm estimates that the relative error in the sum of
             squares is at most TOL.

    INFO = 2 Algorithm estimates that the relative error between X and
             the solution is at most TOL.

    INFO = 3 Conditions for INFO = 1 and INFO = 2 both hold.

    INFO = 4 FVEC is orthogonal to the columns of the Jacobian to
             machine precision.

    INFO = 5 Number of calls to FCN has reached or exceeded 200*(N+1).

    INFO = 6 TOL is too small.  No further reduction in the sum of
             squares is possible.

    INFO = 7 TOL is too small.  No further improvement in the
             approximate solution X is possible.
  */

  return from_double1 (n, x);
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

#ifdef WITH_FFTW_THREADS
  /* FFTW threads: */
  if (nthreads)
    fftw_cleanup_threads ();
#endif
}



/*
  Calling this will define a  few bgy3d-*, hnc3d-*, vec-*, and state-*
  gsubrs introduced above.  This  callback is run by Guile interpreter
  at the latest when the module is imported/compiled.  See the call to
  scm_c_define_module() below.
*/
static void module_init (void* unused)
{
  (void) unused;

  /*
    If Scheme executes  this code inside a module  (which we do), then
    all these gsubrs  will be module procedures available  only in the
    module itself  or by an  explicit (use-modules ...). To  make them
    usable outside of the module one needs to export them. EXPORT() is
    a macro that does both.
  */
  EXPORT ("hnc3d-run-solvent/c", 1, 0, 0, guile_hnc3d_solvent);
  EXPORT ("hnc3d-run-solute/c", 3, 0, 0, guile_hnc3d_solute);
  EXPORT ("bgy3d-run-solvent/c", 1, 0, 0, guile_bgy3d_solvent);
  EXPORT ("bgy3d-run-solute/c", 3, 0, 0, guile_bgy3d_solute);
  EXPORT ("bgy3d-pot-interp", 2, 0, 0, guile_pot_interp);
  EXPORT ("bgy3d-pot-destroy", 1, 0, 0, guile_pot_destroy);
  EXPORT ("bgy3d-restart-destroy", 1, 0, 0, guile_restart_destroy);
  EXPORT ("comm-rank", 0, 0, 0, guile_comm_rank);
  EXPORT ("comm-size", 0, 0, 0, guile_comm_size);
  EXPORT ("comm-set-parallel!", 1, 0, 0, guile_comm_set_parallel_x);
#if SCM_MAJOR_VERSION > 1
  EXPORT ("comm-bcast!", 4, 0, 0, guile_comm_bcast_x);
#endif
  EXPORT ("rism-solvent/c", 1, 0, 0, guile_rism_solvent);
  EXPORT ("rism-solute/c", 2, 1, 0, guile_rism_solute);
  EXPORT ("rism-self-energy/c", 2, 0, 0, guile_rism_self_energy);
  EXPORT ("least-squares", 2, 0, 0, guile_least_squares);
  EXPORT ("bgy3d-test", 3, 0, 0, guile_test);

  /* Define SMOBs: */
  guile_init_state_type ();
  guile_init_vec_type ();

  /* Export force-field primitive shapes used in C: */
  bgy3d_force_field_init ();
}


/* Public API: */
void
bgy3d_guile_init (int argc, char **argv)
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
  assert (sizeof (void (*)()) == sizeof (void*));

  /*
    PETSC  will call  MPI_Init() if  it has  not been  called already.
    PetscFinalize() will  call MPI_Finalize() in that  case also.  MPI
    may choose  to rewrite the command  line, do it  early. Petsc does
    not rewrite argv.  Guile will not understand Petsc flags.
  */
  PetscInitialize (&argc, &argv, NULL, helptext);

  /* Good to know, not required: */
  assert (PETSC_COMM_WORLD == MPI_COMM_WORLD);
  assert (PETSC_COMM_SELF == MPI_COMM_SELF);

  /*
    Declared  in  bgy3d.h,  defined  in bgy3d.c  to  be  MPI_COMM_NULL
    initially.  Here it  is set  for  the first  time. Determines  the
    initial parallelization mode.
  */
  comm_world_petsc = PETSC_COMM_WORLD;

#ifdef WITH_FFTW_THREADS
  if (nthreads)
    {
      /* FFTW threads: */
      const int ok = fftw_init_threads ();
      assert (ok);

      /* FIXME: with MPI parallelization this should be 1, I think: */
      fftw_plan_with_nthreads (nthreads);
    }
#endif

  /* Add  an  exit handler  that  calls  PetscFinalize(). Executed  by
     exit() according to POSIX: */
  atexit (finalize);

  /* Make Petsc abort when it encounters an error: */
  PetscPushErrorHandler (PetscAbortErrorHandler, NULL);

  /*
    FIXME: Set global verbosity early  enough.  This is the only short
    option!   Use   long-options  prefixed   by  "--"  as   the  usual
    convention.
  */
  verbosity = 0; // bgy3d_getopt_test ("-v"); /* extern */

  /* Numeric value, if supplied, will overwrite the on/off "-v"
     flag: */
  // bgy3d_getopt_int ("verbosity", &verbosity);

  /*
   Note that  the names that would  be defined here were  put into the
   private name space of (guile-user) module. For example if you do

     scm_c_define_gsubr ("module-init", 0, 0, 0, module_init);

   (after  adapting interface  of module_init()  accordingly)  then in
   order to  call that  from scheme  you may need  to "steal"  it from
   there by dereferencing as in e.g.

     (@@ (guile-user) module-init).

   Instead we define an internal module (guile bgy3d internal) and put
   all definitions here. When  initialization is actually performed is
   left  to the  interpreter.  In  this way  there  no warnings  about
   "possibly undefined symbols" at compilation stage with Guile 2.0.
  */

  /* scm_c_define_module (const char *name, void (*init)(void *), void
     *data) */
  scm_c_define_module ("guile bgy3d internal", module_init, NULL);
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

/* FIXME: Prototype in bgy3d.h */
void misc_error (const char *loc, const char *msg) /* noreturn */
{
  scm_misc_error (loc, msg, SCM_EOL); /* longjmp! */
}

/* Public API: */
void
bgy3d_molmech (int n, double x[n][3], double *e, double g[n][3])
{
  SCM fg = scm_fluid_ref (lookup ("guile bgy3d", "*server*"));

  /* FIXME: this "if" is questionable, should we rather fail? */
  if (scm_is_true (fg))
    {
      scm_display (fg, scm_current_output_port());
      scm_newline (scm_current_output_port());

      SCM eg = scm_call_1 (fg, from_double2 (n, 3, x));

      scm_display (eg, scm_current_output_port());
      scm_newline (scm_current_output_port());
      /* assert (scm_c_nvalues (eg) == 2); */

      SCM ex = scm_c_value_ref (eg, 0);
      SCM gx = scm_c_value_ref (eg, 1);

      *e = scm_to_double (ex);
      to_double2 (gx, n, 3, g);
    }
  else
    {
      *e = 0.0;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < 3; j++)
          g[i][j] = 0.0;
    }
}
