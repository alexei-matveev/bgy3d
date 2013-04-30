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
void bgy3d_vec_aliases_destroy1 (Vec X, int m, Vec x[m]);

void bgy3d_vec_aliases_create2 (Vec X, int m, Vec x[m][m]);
void bgy3d_vec_aliases_destroy2 (Vec X, int m, Vec x[m][m]);

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


static inline Vec vec_pop (DA da)
{
  Vec work;
  DAGetGlobalVector (da, &work);

  return work;
}


static inline void vec_push (DA da, Vec *work)
{
  DARestoreGlobalVector (da, work);
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
  DAGetInfo (da, &dim, &N[0], &N[1], &N[2],
             NULL, NULL, NULL, NULL, NULL,
             NULL,
#if PETSC_VERSION >= 30200
             NULL, NULL,
#endif
             NULL);
  assert (dim == 3);
}


/* Degrees of freedom, usually 1. For complex vectors -- 2: */
static inline int da_dof (const DA da)
{
  int n;
  DAGetInfo (da, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
             &n, NULL,
             NULL,
#if PETSC_VERSION >= 30200
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


/* Dont forget to vec_restore_array()! */
static inline real* vec_get_array (Vec x)
{
  real *x_;
  VecGetArray (x, &x_);
  return x_;
}

/* Emulates  the behaviour  of the  VecRestoreArray() in  recent PETSC
   versions. */
static inline void vec_restore_array (Vec x, real **x_)
{
  VecRestoreArray (x, x_);
  /* FIXME: More recent PETSC versions nullify the pointer by
     themselves. */
  *x_ = NULL;
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
static inline void bgy3d_vec_map2 (Vec zs, real (*f)(real x, real y), Vec xs, Vec ys)
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


/* ys  = map (f,  xs).  Should  also work  with aliased  arguments for
   in-place transform: */
static inline void bgy3d_vec_fft_map1 (Vec y, /* out */
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
static inline void bgy3d_vec_fft_map2 (Vec z, /* out */
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
static inline void bgy3d_vec_fft_map3 (Vec w, /* out */
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
  DACreate3d (PETSC_COMM_WORLD,
#if PETSC_VERSION >= 30200
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
