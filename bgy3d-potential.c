/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/*
  To use this interface for iterating over the field Vec vec:

  1. Create the context pointer:

     Context *iter = bgy3d_pot_create (BHD, vec);

  2. Iterate  over the grid  doing something usefull with  the values,
     eg. computing the integral:

       while (bgy3d_pot_get_value (iter, chunk_size, x, v, &n))
         for (int i = 0; i < n; i++)
           ... do something with coordinates x[i] and values v[i] ...

     Note that the values from Vec  vec in v[] are sacaled by h^3, the
     grid weight.

  3. You  can  evaluate  potential values  v[n]  at arbitrary  points,
     x[n][3].  Beware of periodic  wrap around and interpolation costs
     of O(N^3) per point:

       bgy3d_pot_interp (iter, n, x, v)

  4. Remember to clean the memory:

     bgy3d_pot_destroy (iter);

  Vec vec does not need to be  a potential, it can be any field on the
  grid.
*/

#include "bgy3d.h"
#include "bgy3d-vec.h"          /* vec_ref(), ... */
#include "bgy3d-fftw.h"         /* bgy3d_fft_interp() */
#include "bgy3d-potential.h"

/*
  Context to mask the explicit form of PETSC stuff.  The corresponding
  typedef struct Context Context is in the header file.
*/
struct Context {
  DA da, dc;                 /* real and complex array descriptors */
  Mat fft_mat;               /* FFT matrix for interpolation */
  Vec v, v_fft;              /* refs to real and complex vectors */
  PetscScalar ***v_;    /* v_[k][j][i] points to the real data */
  int ijk;              /* linarized index for local (k, j, i) */
  real h[3];            /* mesh sizes */
  real interval[2];     /* ? */
  int i0, j0, k0;       /* corner of local grid */
  int ni, nj, nk;       /* local grid shape */
};

/* n = N * j + i, return i and j */
static void divmod (const int n, const int N, int *j, int *i)
{
  *j = n / N;
  *i = n % N;
}


/*
  Put the memory  allocation here.  This has be to  called from C side
  since  we don't  want to  allocate  memory for  vector from  fortran
  return the pointer to context after memory allocation.

  Vec v is saved in the returned iterator context by reference, not by
  copying.  User  code should  not change the  contents of  the vector
  while  iterating over  it.   It may  though bgy3d_vec_destroy()  the
  vector right  after calling this  function as we will  increment the
  refcount.
*/
Context* bgy3d_pot_create (const State *BHD, Vec v)
{
  /* memory for context pointer */
  Context *s = malloc (sizeof *s);

  /* Increment refcounts  and save a  ref in the  Context. DADestroy()
     them! */
  s->da = da_ref (BHD->da);
  s->dc = da_ref (BHD->dc);

  /* This     will     only      be     used     for     trigonometric
     interpolation. MatDestroy() it! */
  s->fft_mat = mat_ref (BHD->fft_mat);

  /*
    Do not copy the input vector, save a reference instead, but also
    Uincrement the refcount. We will have to bgy3d_vec_destroy() it
    too:
  */
  s->v = vec_ref (v);

  /* The first  time interpolation is requested we  put here something
     more usefull: */
  s->v_fft = NULL;

  /* From  now  on  v_[k][j][i]  can   be  used  to  access  a  vector
     element: */
  DAVecGetArray (s->da, s->v, &s->v_);

  /* Get local portion of the grid */
  DAGetCorners (s->da, &s->i0, &s->j0, &s->k0, &s->ni, &s->nj, &s->nk);

  /* copy N3 and h3 from PD */
  FOR_DIM
    s->h[dim] = BHD->PD->h[dim];

  /* copy interval */
  s->interval[0] = BHD->PD->interval[0];
  s->interval[1] = BHD->PD->interval[1];

  /* Initalize counter: */
  s->ijk = 0;

  return s;
}

void bgy3d_pot_interp (Context *s, int n, /* const */ real x[n][3], real v[n])
{
  /* Prepare Fourier coefficients, if not has been already done: */
  if (s->v_fft == NULL)
    {
      s->v_fft = bgy3d_vec_create (s->dc); /* complex */
      MatMult (s->fft_mat, s->v, s->v_fft);
    }

  const real off = s->interval[0];

  /* Translate site coordinates into real grid coordinates where
     integer values correspond to the grid nodes: */
  real y[n][3];
  for (int i = 0; i < n; i++)
    FOR_DIM
      y[i][dim] = (x[i][dim] - off) / s->h[dim];

  /* Trigonometric interpolation: */
  bgy3d_fft_interp (s->fft_mat, s->v_fft, n, y, v);
}


/* Value fetch interface. Returns  number of actually delivered points
   in *np <= n. Returns false on termination. */
bool bgy3d_pot_get_value (Context *s, int n, real x[n][3], real v[n], int *np)
{
  assert (n != 0);            /* need to decide how to handle that! */

  /* How many elements we have: */
  const int local_size = s->ni * s->nj * s->nk;
  const real dV = s->h[0] * s->h[1] * s->h[2];

  /* generating coordinates */
  int p = 0;
  while (p < n && s->ijk < local_size)
    {
      /*
        Get coordinatens i, j, and k within the local block:

        ijk = k * NI * NJ + ij
        ij = j * NI + i
      */
      int ij, i, j, k;
      divmod (s->ijk, s->ni * s->nj, &k, &ij);
      divmod (ij, s->ni, &j, &i);

      /* Add corner coordinates: */
      i += s->i0;
      j += s->j0;
      k += s->k0;

      x[p][0] = i * s->h[0] + s->interval[0]; /* x */
      x[p][1] = j * s->h[1] + s->interval[0]; /* y */
      x[p][2] = k * s->h[2] + s->interval[0]; /* z */
      v[p] = s->v_[k][j][i] * dV;             /* value * weight */

      /* update counters */
      s->ijk++;
      p++;
    }

  /* Reset  counter  to  original   point  once  we  fetched  all  the
     values. FIXME: What if the user supplies n == 0? */
  if (p == 0)
    s->ijk = 0;

  /* Also  return   the  number  of  points  and   a  flag  indicating
     termination: */
  *np = p;

  return p > 0;
}

/* clean up the memory for public *pcontext */
void bgy3d_pot_destroy (Context *s)
{
  /* Required complement of DAVecGetArray(): */
  DAVecRestoreArray (s->da, s->v, &s->v_);

  /* Decrement  refcounts and  destroy  them if  the refcount  reached
     zero: */
  DADestroy (s->da);
  DADestroy (s->dc);
  MatDestroy (&s->fft_mat);
  bgy3d_vec_destroy (&s->v);

  /* Only if interpolation was really used: */
  if (s->v_fft)
    bgy3d_vec_destroy (&s->v_fft);

  /* free the whole context */
  free (s);
}

/* Example usage: */
void bgy3d_pot_test (const State *BHD, Vec vec)
{
  PetscPrintf (PETSC_COMM_WORLD, "Test to potential interface\n");

  /* Make an iterator: */
  Context *s = bgy3d_pot_create (BHD, vec);

  const int chunk_size = 120;
  real v[chunk_size];
  real x[chunk_size][3];
  int nact;

  /* Iterate several times to test wrap-around: */
  for (int times = 0; times < 3; times++)
    {
      /* Calculate moments for tests. First initializing local sums: */
      real m0 = 0.0;
      real m1[3] = {0.0, 0.0, 0.0};
      while (bgy3d_pot_get_value (s, chunk_size, x, v, &nact))
        for (int i = 0; i < nact; i++)
          {
            m0 += v[i];
            for (int j = 0; j < 3; j++)
              m1[j] += x[i][j] * v[i];
          }

      /* Broadcast results of each worker to total sums: */
      bgy3d_comm_allreduce (1, &m0);
      bgy3d_comm_allreduce (3, m1);

      const real L = BHD->PD->interval[1] - BHD->PD->interval[0];
      const real V = L * L * L;
      PetscPrintf (PETSC_COMM_WORLD, "Moments divided by cell volume V = %lf: \n", V);
      PetscPrintf (PETSC_COMM_WORLD, "<1 * v> = %lf\n", m0 / V);
      PetscPrintf (PETSC_COMM_WORLD, "<x * v> = %lf\n", m1[0] / V);
      PetscPrintf (PETSC_COMM_WORLD, "<y * v> = %lf\n", m1[1] / V);
      PetscPrintf (PETSC_COMM_WORLD, "<z * v> = %lf\n", m1[2] / V);
    }

  bgy3d_pot_destroy (s);
}
