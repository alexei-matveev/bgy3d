/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/* FIXME: any better way? */
#include <complex.h>

void vec_rtab (const State *HD, int n, const real rtab[n], real dr,
               Vec v);          /* out */

void vec_ktab (const State *HD, int n, const real ktab[n], real dk,
               Vec v_fft);      /* out */

real bgy3d_vec_mix (Vec dg, Vec dg_new, real a, Vec work);

void bgy3d_vec_save (const char file[], const Vec vec);
void bgy3d_vec_save1 (const char *format, int m, const Vec g[m]);
void bgy3d_vec_save2 (const char *format, int m, /* const */ Vec g2[m][m]);

Vec bgy3d_vec_load (const char file[]); /* Creates a new Vec */

void bgy3d_vec_read (const char file[], Vec vec);  /* Fills existing Vec */
void bgy3d_vec_read1 (const char *format, int m, const Vec g[m]);
void bgy3d_vec_read2 (const char *format, int m, /* const */ Vec g2[m][m]);

void bgy3d_vec_save_ascii (const char file[], const Vec vec);
void bgy3d_vec_save_ascii1 (const char *format, int m, const Vec vec[m]);
void bgy3d_vec_save_ascii2 (const char *format, int m, /* const */ Vec vec[m][m]);


void bgy3d_vec_read_radial (const DA da, const ProblemData *PD,
                            const char *file, Vec g2);
void bgy3d_vec_read_radial2 (const DA da, const ProblemData *PD,
                             const char *format, int m, /* const */ Vec g2[m][m]);

void bgy3d_moments (const State *BHD, Vec v, real q[1], real d[3], real Q[3][3]);

void bgy3d_vec_fft_trans (const DA dc, const int N[static 3], Vec v);


static inline
Vec vec_duplicate (const Vec x)
{
  Vec y;
  VecDuplicate (x, &y);
  return y;
}


static inline
Vec vec_create (const DA da)
{
  Vec x;
  DMCreateGlobalVector (da, &x);
  return x;
}


/* Petsc  3.2 changed the  interface of  XXXDestroy() methods  so that
   they take the pointer to a Petsc object and nullify it: */
static inline
void vec_destroy (Vec *g)
{
  /* Since Petsc 3.2, VecDestroy() also zeroed out the buffer and
   * takes the address instead of vector as input argument */
  VecDestroy (g);
  /* FIXME: only needed before Petsc 3.2 */
  *g = NULL;
}


static inline
void vec_create1 (const DA da, int m, Vec g[m])
{
  for (int i = 0; i < m; i++)
    g[i] = vec_create (da);
}


static inline
void vec_destroy1 (int m, Vec g[m])
{
  for (int i = 0; i < m; i++)
    vec_destroy (&g[i]);
}


/* Allocates g[m][m] with g[j][i] being aliased to g[i][j]: */
static inline
void vec_create2 (const DA da, int m, Vec g[m][m])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      g[j][i] = g[i][j] = vec_create (da);
}


static inline
void vec_destroy2 (int m, Vec g[m][m])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        assert (g[i][j] == g[j][i]);
        vec_destroy (&g[i][j]);
        g[j][i] = NULL;
      }
}


static inline Vec vec_pop (DA da)
{
  Vec work;
  DMGetGlobalVector (da, &work);

  return work;
}


static inline void vec_push (DA da, Vec *work)
{
  DMRestoreGlobalVector (da, work);
  *work = NULL;
}


/* Increment  the  reference count,  and  return  the descriptor,  the
   caller is obliged to DADestroy() it. Do not abuse!. */
static inline DA da_ref (DA da)
{
  PetscObjectReference ((PetscObject) da);
  return da;
}

/* Increment the reference count, and return the vector, the caller is
   obliged to VecDestroy() it. Do not abuse!. */
static inline Vec vec_ref (Vec x)
{
  PetscObjectReference ((PetscObject) x);
  return x;
}

/* Increment the reference count, and return the matrix, the caller is
   obliged to MatDestroy() it. Do not abuse!. */
static inline Mat mat_ref (Mat m)
{
  PetscObjectReference ((PetscObject) m);
  return m;
}


/* Shape of the grid: */
static inline void da_shape (const DA da, int N[static 3])
{
  int dim;
  DMDAGetInfo (da, &dim, &N[0], &N[1], &N[2],
             NULL, NULL, NULL, NULL, NULL,
             NULL,
#if PETSC_VERSION >= VERSION(3, 2)
             NULL, NULL,
#endif
             NULL);
  assert (dim == 3);
}


/* Degrees of freedom, usually 1. For complex vectors -- 2: */
static inline int da_dof (const DA da)
{
  int n;
  DMDAGetInfo (da, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
             &n, NULL,
             NULL,
#if PETSC_VERSION >= VERSION(3, 2)
             NULL, NULL,
#endif
             NULL);
  return n;
}


/*
  To create  a new  Vec one needs  the local  size. The total  size is
  computable.   This is  how  to get  the  local size  from the  array
  descriptor:
*/
static inline int da_local_size (const DA da)
{
  /* Get dimensions and other vector properties: */
  int x[3], n[3];
  DMDAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  const int dof = da_dof (da);  /* dof == 2 for complex Vecs */

  return dof * n[2] * n[1] * n[0];
}


static inline Vec vec_from_array (int n, real x_[n])
{
  Vec x;
#if PETSC_VERSION >= VERSION(3, 3)
  /*
    Interface changed  in 3.3!  The  second argument, PetscInt  bs, is
    the block size that  is apparently used for VecSetValuesBlocked().
    FIXME: literal constant here!
  */
  VecCreateMPIWithArray (comm_world_petsc, 1, n, PETSC_DECIDE, x_, &x);
#else
  VecCreateMPIWithArray (comm_world_petsc, n, PETSC_DECIDE, x_, &x);
#endif
  return x;
}


/* Maybe Vec -> Maybe Ptr.  Dont forget to vec_restore_array()! */
static inline real* vec_get_array (Vec x)
{
  /* Some Vec y[m][m] are NULL on the diagonal. This hack is to handle
     them transparently: */
  if (x == NULL)
    return NULL;

  real *x_;
  VecGetArray (x, &x_);
  return x_;
}


/*
  Maybe Vec -> Ptr (Maybe Ptr) -> Void.  Nullifies *x_ which otherwise
  would be  dangling semantically.  This behaviour  is consistent with
  that of the VecRestoreArray() in recent PETSC versions.
*/
static inline void vec_restore_array (Vec x, real **x_)
{
  /* The couterpart to the hack in vec_get_array(): */
  if (x == NULL && *x_ == NULL)
    return;

  VecRestoreArray (x, x_);
  /* FIXME:  More  recent  PETSC   versions  nullify  the  pointer  by
     themselves. */
  *x_ = NULL;
}


static inline void vec_get_array2 (int m, Vec x[m][m], real *x_[m][m])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        x_[i][j] = x_[j][i] = vec_get_array (x[i][j]);
        assert (x[i][j] == x[j][i]);
      }
}


static inline void vec_restore_array2 (int m, Vec x[m][m], real *x_[m][m])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        assert (x[i][j] == x[j][i]);
        assert (x_[i][j] == x_[j][i]);
        vec_restore_array (x[i][j], &x_[i][j]);
        x_[j][i] = NULL;
      }
}


static inline int vec_local_size (Vec x)
{
  int n;
  VecGetLocalSize (x, &n);
  return n;
}

static inline int vec_size (Vec x)
{
  int n;
  VecGetSize (x, &n);
  return n;
}

static inline real vec_sum (Vec x)
{
  real sum;
  VecSum (x, &sum);
  return sum;
}

static inline real vec_avg (Vec x)
{
  return vec_sum (x) / vec_size (x);
}

static inline real vec_norm (Vec x)
{
  real norm;
  VecNorm (x, NORM_INFINITY, &norm);
  return norm;
}

static inline real vec_dot (Vec x, Vec y)
{
  real dot;
  VecDot (x, y, &dot);
  return dot;
}

static inline real vec_max (Vec x)
{
  real max;
  VecMax (x, NULL, &max);      /* dont want location */
  return max;
}

static inline real vec_min (Vec x)
{
  real min;
  VecMin (x, NULL, &min);      /* dont want location */
  return min;
}


/*
  The  following vec_map1(), vec_map2(),  vec_map3() should  be better
  inlined as they use a scalar function f().

    ys = map (f, xs)

  Should also work with aliased arguments for in-place transform.
*/
static inline void
vec_map1 (Vec ys, real (*f)(real x), Vec xs)
{
  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);

  local real *xs_ = vec_get_array (xs);
  local real *ys_ = vec_get_array (ys);

  for (int i = 0; i < n; i++)
    ys_[i] = f (xs_[i]);

  vec_restore_array (xs, &xs_);
  vec_restore_array (ys, &ys_);
}

/* zs = map (f, xs, ys).   Should also work with aliased arguments for
   in-place transform: */
static inline void
vec_map2 (Vec zs, real (*f)(real x, real y), Vec xs, Vec ys)
{
  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);
  assert (vec_local_size (zs) == n);

  local real *xs_ = vec_get_array (xs);
  local real *ys_ = vec_get_array (ys);
  local real *zs_ = vec_get_array (zs);

  for (int i = 0; i < n; i++)
    zs_[i] = f (xs_[i], ys_[i]);

  vec_restore_array (xs, &xs_);
  vec_restore_array (ys, &ys_);
  vec_restore_array (zs, &zs_);
}


/* ws = map (f, xs, ys,  zs).  Should also work with aliased arguments
   for in-place transform. Consider using vec_app4() instead. */
deprecated static inline void
vec_map3 (Vec ws, real (*f)(real x, real y, real z),
          Vec xs, Vec ys, Vec zs)
{
  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);
  assert (vec_local_size (zs) == n);
  assert (vec_local_size (ws) == n);

  local real *xs_ = vec_get_array (xs);
  local real *ys_ = vec_get_array (ys);
  local real *zs_ = vec_get_array (zs);
  local real *ws_ = vec_get_array (ws);

  for (int i = 0; i < n; i++)
    ws_[i] = f (xs_[i], ys_[i], zs_[i]);

  vec_restore_array (xs, &xs_);
  vec_restore_array (ys, &ys_);
  vec_restore_array (zs, &zs_);
  vec_restore_array (ws, &ws_);
}


/*
  These vec_app?() functions use arrays  and do not need to be inlined
  solely  for performance reasons.   But I  hate prefixing  them. This
  code does not dictate what is in- and what is output. In fact all of
  them may  be modified.  But should  they?  One convention  is to use
  the first argument as output, the rest for input:
*/
static inline void
vec_app3 (void (*f)(int n, real x0[n], real x1[n], real x2[n]),
          Vec x0, Vec x1, Vec x2)
{
  const int n = vec_local_size (x0);
  assert (vec_local_size (x1) == n);
  assert (vec_local_size (x2) == n);

  local real *x0_ = vec_get_array (x0);
  local real *x1_ = vec_get_array (x1);
  local real *x2_ = vec_get_array (x2);

  f (n, x0_, x1_, x2_);

  vec_restore_array (x0, &x0_);
  vec_restore_array (x1, &x1_);
  vec_restore_array (x2, &x2_);
}


static inline void
vec_app4 (void (*f)(int n, real x0[n], real x1[n], real x2[n], real x3[n]),
          Vec x0, Vec x1, Vec x2, Vec x3)
{
  const int n = vec_local_size (x0);
  assert (vec_local_size (x1) == n);
  assert (vec_local_size (x2) == n);
  assert (vec_local_size (x3) == n);

  local real *x0_ = vec_get_array (x0);
  local real *x1_ = vec_get_array (x1);
  local real *x2_ = vec_get_array (x2);
  local real *x3_ = vec_get_array (x3);

  f (n, x0_, x1_, x2_, x3_);

  vec_restore_array (x0, &x0_);
  vec_restore_array (x1, &x1_);
  vec_restore_array (x2, &x2_);
  vec_restore_array (x3, &x3_);
}


static inline void
vec_app5 (void (*f)(int n, real x0[n], real x1[n], real x2[n], real x3[n], real x4[n]),
          Vec x0, Vec x1, Vec x2, Vec x3, Vec x4)
{
  const int n = vec_local_size (x0);
  assert (vec_local_size (x1) == n);
  assert (vec_local_size (x2) == n);
  assert (vec_local_size (x3) == n);
  assert (vec_local_size (x4) == n);

  local real *x0_ = vec_get_array (x0);
  local real *x1_ = vec_get_array (x1);
  local real *x2_ = vec_get_array (x2);
  local real *x3_ = vec_get_array (x3);
  local real *x4_ = vec_get_array (x4);

  f (n, x0_, x1_, x2_, x3_, x4_);

  vec_restore_array (x0, &x0_);
  vec_restore_array (x1, &x1_);
  vec_restore_array (x2, &x2_);
  vec_restore_array (x3, &x3_);
  vec_restore_array (x4, &x4_);
}


/* FIXME: get rid of chem. pot. form as a functional of many arguments
   x, h, c, cl, ... And you wont need this monster: */
static inline void
vec_app8 (void (*f)(int n,
                    real x0[n], real x1[n], real x2[n], real x3[n],
                    real x4[n], real x5[n], real x6[n], real x7[n]),
          Vec x0, Vec x1, Vec x2, Vec x3,
          Vec x4, Vec x5, Vec x6, Vec x7)
{
  const int n = vec_local_size (x0);
  assert (vec_local_size (x1) == n);
  assert (vec_local_size (x2) == n);
  assert (vec_local_size (x3) == n);
  assert (vec_local_size (x4) == n);
  assert (vec_local_size (x5) == n);
  assert (vec_local_size (x6) == n);
  assert (vec_local_size (x7) == n);

  local real *x0_ = vec_get_array (x0);
  local real *x1_ = vec_get_array (x1);
  local real *x2_ = vec_get_array (x2);
  local real *x3_ = vec_get_array (x3);
  local real *x4_ = vec_get_array (x4);
  local real *x5_ = vec_get_array (x5);
  local real *x6_ = vec_get_array (x6);
  local real *x7_ = vec_get_array (x7);

  f (n, x0_, x1_, x2_, x3_, x4_, x5_, x6_, x7_);

  vec_restore_array (x0, &x0_);
  vec_restore_array (x1, &x1_);
  vec_restore_array (x2, &x2_);
  vec_restore_array (x3, &x3_);
  vec_restore_array (x4, &x4_);
  vec_restore_array (x5, &x5_);
  vec_restore_array (x6, &x6_);
  vec_restore_array (x7, &x7_);
}


/* ys  = map (f,  xs).  Should  also work  with aliased  arguments for
   in-place transform: */
static inline void vec_fft_map1 (Vec y, /* out */
                                 complex (*f)(complex x),
                                 Vec x) /* in */
{
  const int n = vec_local_size (x);
  assert (vec_local_size (y) == n);
  assert (n % 2 == 0);

  local complex *x_ = (complex*) vec_get_array (x);
  local complex *y_ = (complex*) vec_get_array (y);

  for (int i = 0; i < n / 2; i++)
    y_[i] = f (x_[i]);

  vec_restore_array (x, (void*) &x_);
  vec_restore_array (y, (void*) &y_);
}


/* zs = map (f, xs, ys).   Should also work with aliased arguments for
   in-place transform: */
static inline void vec_fft_map2 (Vec z, /* out */
                                 complex (*f)(complex x, complex y),
                                 Vec x, Vec y) /* in */
{
  const int n = vec_local_size (x);
  assert (vec_local_size (y) == n);
  assert (vec_local_size (z) == n);
  assert (n % 2 == 0);

  local complex *x_ = (complex*) vec_get_array (x);
  local complex *y_ = (complex*) vec_get_array (y);
  local complex *z_ = (complex*) vec_get_array (z);

  for (int i = 0; i < n / 2; i++)
    z_[i] = f (x_[i], y_[i]);

  vec_restore_array (x, (void*) &x_);
  vec_restore_array (y, (void*) &y_);
  vec_restore_array (z, (void*) &z_);
}


/* ws = map (f, xs, ys,  zs).  Should also work with aliased arguments
   for in-place transform: */
static inline void vec_fft_map3 (Vec w, /* out */
                                 complex (*f)(complex x, complex y, complex z),
                                 Vec x, Vec y, Vec z) /* in */
{
  const int n = vec_local_size (x);
  assert (vec_local_size (y) == n);
  assert (vec_local_size (z) == n);
  assert (vec_local_size (w) == n);
  assert (n % 2 == 0);

  local complex *x_ = (complex*) vec_get_array (x);
  local complex *y_ = (complex*) vec_get_array (y);
  local complex *z_ = (complex*) vec_get_array (z);
  local complex *w_ = (complex*) vec_get_array (w);

  for (int i = 0; i < n / 2; i++)
    w_[i] = f (x_[i], y_[i], z_[i]);

  vec_restore_array (x, (void*) &x_);
  vec_restore_array (y, (void*) &y_);
  vec_restore_array (z, (void*) &z_);
  vec_restore_array (w, (void*) &w_);
}


/*
  Tabulate v = f(r) with origin at the grid center. Here r[3] are the
  coordinates of the grid point (small r are in the grid center).
*/
static inline
void vec_rmap3 (const State *BHD, real (*f)(const real r[3]), Vec v)
{
  const real *L = BHD->PD->L;   /* [3] */
  const real *h = BHD->PD->h;   /* [3] */

  real ***v_;
  DMDAVecGetArray (BHD->da, v, &v_);

  int n[3], a[3], i[3];
  DMDAGetCorners (BHD->da, &a[0], &a[1], &a[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = a[2]; i[2] < a[2] + n[2]; i[2]++)
    for (i[1] = a[1]; i[1] < a[1] + n[1]; i[1]++)
      for (i[0] = a[0]; i[0] < a[0] + n[0]; i[0]++)
        {
          real r[3];
          FOR_DIM
            r[dim] = i[dim] * h[dim] - L[dim] / 2;

          v_[i[2]][i[1]][i[0]] = f (r);
        }
  DMDAVecRestoreArray (BHD->da, v, &v_);
}


/*
  Tabulate v_fft = f(k) with origin  at the grid corner. Here k[3] are
  the  coordinates  of the  k-grid  point (small  k  are  in the  grid
  corner).
*/
static inline
void vec_kmap3 (const State *BHD, complex (*f)(const real k[3]), Vec v_fft)
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* [3] */

  real dk[3];                   /* k-mesh spacing */
  FOR_DIM
    dk[dim] = 2 * M_PI / PD->L[dim];

  /* Get local portion of the k-grid */
  int a[3], n[3], i[3];
  DMDAGetCorners (BHD->dc, &a[0], &a[1], &a[2], &n[0], &n[1], &n[2]);

  complex ***v_fft_;
  DMDAVecGetArray (BHD->dc, v_fft, &v_fft_);

   /* loop over local portion of grid */
  for (i[2] = a[2]; i[2] < a[2] + n[2]; i[2]++)
    for (i[1] = a[1]; i[1] < a[1] + n[1]; i[1]++)
      for (i[0] = a[0]; i[0] < a[0] + n[0]; i[0]++)
        {
          real k[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            k[dim] = KFREQ (i[dim], N[dim]) * dk[dim];

          v_fft_[i[2]][i[1]][i[0]] = f (k);
        }
  DMDAVecRestoreArray (BHD->dc, v_fft, &v_fft_);
}


/* Tabulate v = f(r) with origin at the grid center:  */
static inline
void vec_rmap (const State *BHD, real (*f)(real r), Vec v)
{
  real f3 (const real x[3])
  {
    const real r2 = SQR (x[0]) + SQR (x[1]) + SQR (x[2]);
    return f (sqrt (r2));
  }
  vec_rmap3 (BHD, f3, v);
}


/*
  Tabulate v_fft = f(k) with origin at the grid corner.

  FIXME: we take a square root  here, but in the most interesting case
  the Coulomb depends on k^2  anyway. This potential is also a complex
  number with zero imaginary part.
*/
static inline
void vec_kmap (const State *BHD, complex (*f)(real k), Vec v_fft)
{
  complex f3 (const real k[3])
  {
    const real k2 = SQR (k[0]) + SQR (k[1]) + SQR (k[2]);
    return f (sqrt (k2));
  }
  vec_kmap3 (BHD, f3, v_fft);
}


/* "Integrates" f(v(x), x) with the grid data v(x): */
static inline real vec_integrate (DA da, real (*f)(real v, int i, int j, int k), Vec v)
{
  /* Clear accumulator: */
  real acc = 0.0;

  real ***v_;
  DMDAVecGetArray (da, v, &v_);

  /* Loop over local portion of grid */
  int x[3], n[3], i[3];
  DMDAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        acc += f (v_[i[2]][i[1]][i[0]], i[0], i[1], i[2]);

  DMDAVecRestoreArray (da, v, &v_);

  /* Sum accumulator over workers: */
  comm_allreduce (1, &acc);

  return acc;
}

/* Returns sum (1 - g): */
static inline real vec_hole (Vec g)
{
  /* FIXME: loss of precision possible here: */
  return vec_size (g) - vec_sum (g);
}


static inline
void vec_aliases_create1 (Vec X, int m, Vec x[m])
{
  /*
    The length of Vec X should be divisible by m.  Though in principle
    any Vec  X satisfying this should  be accepted, this  code is only
    used for Vecs created by vec_pack_create1(), see below:
  */
  const int mn = vec_local_size (X);
  const int n = mn / m;
  assert (m * n == mn);

  /*
    This  buf has  enough space  for m  Vecs.  Note  that there  is no
    corresponding vec_restore_array() call in  this scope and the lack
    of the "local"  attribute. The contents of the  long Vec X remains
    "checked  out"  for  the  whole  lifetime of  the  aliases.   This
    lifetime ends  upon call  to vec_aliases_destroy1(). Only  then (a
    copy of) the pointer will be "returned".
  */
  real *buf = vec_get_array (X);

  /* Create m Vecs with the storage from buf: */
  for (int i = 0; i < m; i++)
    x[i] = vec_from_array (n, buf + i * n);
}


static inline
void vec_aliases_create2 (Vec X, int m, Vec x[m][m])
{
  /* The length of Vec X should be divisible by m * (m + 1) / 2!*/
  const int nm2 = vec_local_size (X);
  const int m2 = m * (m + 1) / 2;
  const int n = nm2 / m2;
  assert (n * m2 == nm2);

  /* Enough   space    for   m2    Vecs.   See   also    comments   in
     vec_aliases_create1(): */
  real *buf = vec_get_array (X);

  /* Create m2 Vecs with the storage from buf: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        x[i][j] = x[j][i] = vec_from_array (n, buf);
        buf += n;
      }
}


/*
  This and the next function  should not attempt to free() the storage
  of  the aliases  x[].  It  is owned  by the  longer Vec  X.   We are
  relying on the magic of VecDestroy()  that alone knows how a Vec was
  created --- it  should not free() the storage if  Vec was created by
  vec_from_array().
*/
static inline
void vec_aliases_destroy1 (Vec X, int m, Vec x[m])
{
  vec_destroy1 (m, x);    /* should not free() */

  /*
    The epoch  of accessing the  content of Vec  X via the  aliases is
    over. Signal to PETSC that the  content of the Vec X may have been
    changed,  so  that  it  invalidates eventually  cached  derivative
    values such as the vector norm. It is assumed that vec_get_array()
    is idempotent (returns the same value on succesive calls).
  */
  local real *X_ = vec_get_array (X);
  vec_restore_array (X, &X_);
}


static inline
void vec_aliases_destroy2 (Vec X, int m, Vec x[m][m])
{
  vec_destroy2 (m, x);    /* should not free() */

  /* See comments in vec_aliases_destroy1(): */
  local real *X_ = vec_get_array (X);
  vec_restore_array (X, &X_);
}


/*
  Create   a  vector   m-times  longer   than  the   array  descriptor
  specification. See vec_aliases_create*() for what may happen to such
  Vec later.
*/
static inline
Vec vec_pack_create1 (const DA da, int m)
{
  /* Allocate space for m Vecs: */
  const int mn = m * da_local_size (da);

  return vec_from_array (mn, malloc (mn * sizeof (real)));
}


/* Nearly the same as vec_pack_create1(): */
static inline
Vec vec_pack_create2 (const DA da, int m)
{
  return vec_pack_create1 (da, m * (m + 1) / 2);
}


/* VecDestroy() will  not free the storage  if it was  provided by the
   user. We do it ourselves: */
static inline
void vec_pack_destroy1 (Vec *X)
{
  /* FIXME: should we also vec_restore_array()? */
  free (vec_get_array (*X));    /* free() the whole */

  /* take the address of vector as input since Petsc 3.2 */
  VecDestroy (X);
  /* FIXME: only needed before Petsc 3.2 */
  *X = NULL;
}


/* Same as vec_pack_destroy1() */
static inline
void vec_pack_destroy2 (Vec *X)
{
  vec_pack_destroy1 (X);
}


/*
  Create  descriptor of  a distributed  M x  N x  P array  with lx[m],
  ly[n], and lz[p] being partitions of M, N, and P, respectively.
*/
static inline DA da_create (int dof,
                            int m, /* const */ int lx[m],
                            int n, /* const */ int ly[n],
                            int p, /* const */ int lz[p])
{
  const PetscInt stencil_width = 1;

  /* No need to sum over workers, everyone provides the same input: */
  const int M = isum (m, lx);
  const int N = isum (n, ly);
  const int P = isum (p, lz);

  /*
    As of  Petsc 3 the library  refuses (loudly) to  handle zero range
    for any worker.  This limits the  maximum number of workers by P =
    sum(lz).   The boundary type  was non-periodic  in the  very first
    version. FIXME:  I assume it  is only necessary for  the Dirichlet
    boundary  conditions that run  finite stencil  over the  grid. See
    bgy3d.h for BOUNDARY_TYPE and STENCIL_TYPE macros.
  */
  DA da;
  DMDACreate3d (comm_world_petsc,
#if PETSC_VERSION >= VERSION(3, 2)
              BOUNDARY_TYPE, BOUNDARY_TYPE, BOUNDARY_TYPE,
#else
              BOUNDARY_TYPE,
#endif
              STENCIL_TYPE,
              M, N, P,          /* grid dimensions */
              m, n, p,          /* m = 1, n = 1, p = NP, usually */
              dof,              /* 1 or 2, usually */
              stencil_width,    /* 1 */
              lx, ly, lz,       /* partitions of M, N, and P */
              &da);

  return da;
}
