/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/* FIXME: any better way? */
#include <complex.h>

Vec bgy3d_vec_create (const DA da);
void bgy3d_vec_destroy (Vec g);

void bgy3d_vec_create1 (const DA da, int m, Vec g[m]);
void bgy3d_vec_destroy1 (int m, Vec g[m]);

void bgy3d_vec_create2 (const DA da, int m, Vec g[m][m]);
void bgy3d_vec_destroy2 (int m, Vec g[m][m]);

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


void bgy3d_vec_read_radial2 (const State *BHD,
                             const char *format, int m, /* const */ Vec g2[m][m]);

void bgy3d_vec_moments (const DA da, Vec v,
                        real *q, real *x, real *y, real *z);

/* Increment  the  reference count,  and  return  the descriptor,  the
   caller is obliged to VecDestroy() it. Do not abuse!. */
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
  real *xs_, *ys_;

  VecGetArray (xs, &xs_);
  VecGetArray (ys, &ys_);

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
    real *xs_, *ys_, *zs_;

  VecGetArray (xs, &xs_);
  VecGetArray (ys, &ys_);
  VecGetArray (zs, &zs_);

  const int n = vec_local_size (xs);
  assert (vec_local_size (ys) == n);
  assert (vec_local_size (zs) == n);

  for (int i = 0; i < n; i++)
    zs_[i] = f (xs_[i], ys_[i]);

  VecRestoreArray (xs, &xs_);
  VecRestoreArray (ys, &ys_);
  VecRestoreArray (zs, &zs_);
}

/* zs = map (f, xs, ys).   Should also work with aliased arguments for
   in-place transform: */
static inline void bgy3d_vec_fft_map2 (Vec z, complex (*f)(complex x, complex y), Vec x, Vec y)
{
  real *x_, *y_, *z_;

  VecGetArray (x, &x_);
  VecGetArray (y, &y_);
  VecGetArray (z, &z_);

  const int n = vec_local_size (x);
  assert (vec_local_size (y) == n);
  assert (vec_local_size (z) == n);
  assert (n % 2 == 0);

  complex *xs_ = (complex*) x_;
  complex *ys_ = (complex*) y_;
  complex *zs_ = (complex*) z_;

  for (int i = 0; i < n / 2; i++)
    zs_[i] = f (xs_[i], ys_[i]);

  VecRestoreArray (x, &x_);
  VecRestoreArray (y, &y_);
  VecRestoreArray (z, &z_);
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
