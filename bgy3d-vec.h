/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/* FIXME: any better way? */
#include <complex.h>

Vec bgy3d_vec_duplicate (const Vec x);
Vec bgy3d_vec_create (const DA da);
void bgy3d_vec_destroy (Vec *g);

void bgy3d_vec_create1 (const DA da, int m, Vec g[m]);
void bgy3d_vec_destroy1 (int m, Vec g[m]);

void bgy3d_vec_create2 (const DA da, int m, Vec g[m][m]);
void bgy3d_vec_destroy2 (int m, Vec g[m][m]);

void bgy3d_vec_aliases_create1 (Vec X, int m, Vec x[m]);
void bgy3d_vec_aliases_destroy1 (int m, Vec x[m]);

void bgy3d_vec_aliases_create2 (Vec X, int m, Vec x[m][m]);
void bgy3d_vec_aliases_destroy2 (int m, Vec x[m][m]);

Vec bgy3d_vec_pack_create1 (const DA da, int m);
void bgy3d_vec_pack_destroy1 (Vec *X);

Vec bgy3d_vec_pack_create2 (const DA da, int m);
void bgy3d_vec_pack_destroy2 (Vec *X);

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

void bgy3d_vec_moments (const DA da, Vec v,
                        real *q, real *x, real *y, real *z);

void bgy3d_vec_fft_trans (const DA dc, const int N[static 3], Vec v);


static inline Vec bgy3d_vec_pop (DA da)
{
  Vec work;
  DAGetGlobalVector (da, &work);

  return work;
}


static inline void bgy3d_vec_push (DA da, Vec *work)
{
  DARestoreGlobalVector (da, work);
  *work = NULL;
}


/* Increment  the  reference count,  and  return  the descriptor,  the
   caller is obliged to DADestroy() it. Do not abuse!. */
static inline DA bgy3d_da_ref (DA da)
{
  PetscObjectReference ((PetscObject) da);
  return da;
}

/* Increment the reference count, and return the vector, the caller is
   obliged to VecDestroy() it. Do not abuse!. */
static inline Vec bgy3d_vec_ref (Vec x)
{
  PetscObjectReference ((PetscObject) x);
  return x;
}

/* Increment the reference count, and return the matrix, the caller is
   obliged to MatDestroy() it. Do not abuse!. */
static inline Mat bgy3d_mat_ref (Mat m)
{
  PetscObjectReference ((PetscObject) m);
  return m;
}


/* Shape of the grid: */
static inline void da_shape (const DA da, int N[static 3])
{
  int dim;
  DAGetInfo (da, &dim, &N[0], &N[1], &N[2],
             NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  assert (dim == 3);
}


/* Degrees of freedom, usually 1. For complex vectors -- 2: */
static inline int da_dof (const DA da)
{
  int n;
  DAGetInfo (da, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
             &n, NULL, NULL, NULL);
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
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  const int dof = da_dof (da);  /* dof == 2 for complex Vecs */

  return dof * n[2] * n[1] * n[0];
}


static inline Vec vec_from_array (int n, real x_[n])
{
  Vec x;
  /* FIXME: interface changed in 3.3! */
  VecCreateMPIWithArray (PETSC_COMM_WORLD, n, PETSC_DECIDE, x_, &x);
  return x;
}

static inline real* vec_get_array (Vec x)
{
  real *x_;
  VecGetArray (x, &x_);
  return x_;
}

static inline int vec_local_size (Vec x)
{
  int n;
  VecGetLocalSize (x, &n);
  return n;
}

static inline int bgy3d_vec_size (Vec x)
{
  int n;
  VecGetSize (x, &n);
  return n;
}

static inline real bgy3d_vec_sum (Vec x)
{
  real sum;
  VecSum (x, &sum);
  return sum;
}

static inline real bgy3d_vec_avg (Vec x)
{
  return bgy3d_vec_sum (x) / bgy3d_vec_size (x);
}

static inline real bgy3d_vec_norm (Vec x)
{
  real norm;
  VecNorm (x, NORM_INFINITY, &norm);
  return norm;
}

static inline real bgy3d_vec_dot (Vec x, Vec y)
{
  real dot;
  VecDot (x, y, &dot);
  return dot;
}

static inline real bgy3d_vec_max (Vec x)
{
  real max;
  VecMax (x, NULL, &max);      /* dont want location */
  return max;
}

static inline real bgy3d_vec_min (Vec x)
{
  real min;
  VecMin (x, NULL, &min);      /* dont want location */
  return min;
}


/* ys  = map  (f, xs).  Should also  work with  aliased  arguments for
   in-place transform: */
static inline void bgy3d_vec_map1 (Vec ys, real (*f)(real x), Vec xs)
{
  real *xs_ = vec_get_array (xs);
  real *ys_ = vec_get_array (ys);

  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);

  for (int i = 0; i < n; i++)
    ys_[i] = f (xs_[i]);

  VecRestoreArray (xs, &xs_);
  VecRestoreArray (ys, &ys_);
}

/* zs = map (f, xs, ys).   Should also work with aliased arguments for
   in-place transform: */
static inline void bgy3d_vec_map2 (Vec zs, real (*f)(real x, real y), Vec xs, Vec ys)
{
  real *xs_ = vec_get_array (xs);
  real *ys_ = vec_get_array (ys);
  real *zs_ = vec_get_array (zs);

  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);
  assert (vec_local_size (zs) == n);

  for (int i = 0; i < n; i++)
    zs_[i] = f (xs_[i], ys_[i]);

  VecRestoreArray (xs, &xs_);
  VecRestoreArray (ys, &ys_);
  VecRestoreArray (zs, &zs_);
}


/* ys  = map (f,  xs).  Should  also work  with aliased  arguments for
   in-place transform: */
static inline void bgy3d_vec_fft_map1 (Vec y, /* out */
                                       complex (*f)(complex x),
                                       Vec x) /* in */
{
  complex *x_ = (complex*) vec_get_array (x);
  complex *y_ = (complex*) vec_get_array (y);

  const int n = vec_local_size (x);
  assert (vec_local_size (y) == n);
  assert (n % 2 == 0);

  for (int i = 0; i < n / 2; i++)
    y_[i] = f (x_[i]);

  VecRestoreArray (x, (void*) &x_);
  VecRestoreArray (y, (void*) &y_);
}


/* zs = map (f, xs, ys).   Should also work with aliased arguments for
   in-place transform: */
static inline void bgy3d_vec_fft_map2 (Vec z, /* out */
                                       complex (*f)(complex x, complex y),
                                       Vec x, Vec y) /* in */
{
  complex *x_ = (complex*) vec_get_array (x);
  complex *y_ = (complex*) vec_get_array (y);
  complex *z_ = (complex*) vec_get_array (z);

  const int n = vec_local_size (x);
  assert (vec_local_size (y) == n);
  assert (vec_local_size (z) == n);
  assert (n % 2 == 0);

  for (int i = 0; i < n / 2; i++)
    z_[i] = f (x_[i], y_[i]);

  VecRestoreArray (x, (void*) &x_);
  VecRestoreArray (y, (void*) &y_);
  VecRestoreArray (z, (void*) &z_);
}


/* ws = map (f, xs, ys,  zs).  Should also work with aliased arguments
   for in-place transform: */
static inline void bgy3d_vec_fft_map3 (Vec w, /* out */
                                       complex (*f)(complex x, complex y, complex z),
                                       Vec x, Vec y, Vec z) /* in */
{
  complex *x_ = (complex*) vec_get_array (x);
  complex *y_ = (complex*) vec_get_array (y);
  complex *z_ = (complex*) vec_get_array (z);
  complex *w_ = (complex*) vec_get_array (w);

  const int n = vec_local_size (x);
  assert (vec_local_size (y) == n);
  assert (vec_local_size (z) == n);
  assert (vec_local_size (w) == n);
  assert (n % 2 == 0);

  for (int i = 0; i < n / 2; i++)
    w_[i] = f (x_[i], y_[i], z_[i]);

  VecRestoreArray (x, (void*) &x_);
  VecRestoreArray (y, (void*) &y_);
  VecRestoreArray (z, (void*) &z_);
  VecRestoreArray (w, (void*) &w_);
}


/* "Integrates" f(v(x), x) with the grid data v(x): */
static inline real bgy3d_vec_integrate (DA da, real (*f)(real v, int i, int j, int k), Vec v)
{
  /* Clear accumulator: */
  real acc = 0.0;

  real ***v_;
  DAVecGetArray (da, v, &v_);

  /* Loop over local portion of grid */
  int x[3], n[3], i[3];
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        acc += f (v_[i[2]][i[1]][i[0]], i[0], i[1], i[2]);

  DAVecRestoreArray (da, v, &v_);

  /* Sum accumulator over workers: */
  bgy3d_comm_allreduce (1, &acc);

  return acc;
}

/* Returns sum (1 - g): */
static inline real bgy3d_vec_hole (Vec g)
{
  /* FIXME: loss of precision possible here: */
  return bgy3d_vec_size (g) - bgy3d_vec_sum (g);
}
