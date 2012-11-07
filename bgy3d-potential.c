/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-potential.h"

/* Context to mask the explicit form of PETSC stuff */
typedef struct Context {
  Vec v;    /* vector storing potential values */
  real (*x)[3]; /* (huge)  array for  coordinates. FIXME:  maybe  also use
               vector? */
  int counter;              /* counter indicating how many are left */
  int nmax; /* length  of vector and  array, better  save it  than use
               VecGetSize() from time to time */
} Context;

/* Put the memory allocation here.  This  has be to called from C side
 since we don't want to allocate memory for vector from fortran return the void
 pointer to context after memory allocation */
void* bgy3d_pot_create (DA da, const ProblemData *PD, Vec v)
{
  int i0, j0, k0;
  int ni, nj, nk;

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);
  int m = ni * nj * nk;

  /* memory for context pointer */
  Context *pcontext = malloc(sizeof *pcontext);

  /* allocate memory for coordinate member in context */
  pcontext->x = malloc(m * 3 * sizeof(real));

  /* Get coordinates of the local grid portion: */
  {
    int ijk = 0;
    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          {
            /* coordinates */
            /* pcontext->x[ijk][0] = i * PD->h[0] + PD->interval[0];
               pcontext->x[ijk][1] = j * PD->h[1] + PD->interval[0];
               pcontext->x[ijk][2] = k * PD->h[2] + PD->interval[0]; */

            /* use grid index for test only */
            pcontext->x[ijk][0] = i;
            pcontext->x[ijk][1] = j;
            pcontext->x[ijk][2] = k;
            ijk++;
          }
    assert (ijk == m);
  }

  /* create vector member in context */
  DACreateGlobalVector (da, &pcontext->v);

  /* Copy the value from input vector */
  VecCopy (v, pcontext->v);

  /* set counter number and save vector length */
  pcontext->counter = m;
  pcontext->nmax = m;

  return (void *)pcontext;
}

/* Value fetch interface */
int bgy3d_pot_get_value (void *s, int n, real x[n][3], real v[n])
{
  /* type cast */
  Context *pct = (Context *)s;

  /* head index for context->v and context->x */
  int nstart = pct->nmax - pct->counter;

  /* number of values that would be actually fetched */
  int nact = MIN(n, pct->counter);

  /* fill context->x */
  for (int i = 0; i < nact; i++)
    for ( int j = 0; j < 3; j++)
      x[i][j] = pct->x[nstart + i][j];

  /* array storing indices that would be fetched from context->v */
  int idx[nact];
  for (int i = 0; i < nact; i++)
    idx[i] = nstart + i;

  /* Get values from context->v[nstart : nstart + nact] */
  PetscErrorCode ierr = VecGetValues(pct->v, nact, idx, v);
  assert (!ierr);

  /* update counter */
  pct->counter -= nact;

  /* reset counter to original point once we fetched all the values */
  if (nact == 0)
    pct->counter = pct->nmax;

  return nact;
}

/* clean up the memory for public *pcontext */
void bgy3d_pot_destroy (void* s)
{
  /* type cast, again? */
  Context *pct = (Context *)s;

  /* free the dynamically allocated context->x */
  free(pct->x);

  /* free context->v */
  VecDestroy (pct->v);

  /* free the whole context */
  free(pct);
}

/* Test for interface */
void bgy3d_pot_test (void *s)
{
  const int chunk_size = 120;
  real v[chunk_size];
  real x[chunk_size][3];
  int nact;


  /* calculate moments for tests: */
  real m0 = 0.0;
  real mx = 0.0;
  real my = 0.0;
  real mz = 0.0;
  while ((nact = bgy3d_pot_get_value (s, chunk_size, x, v)))
    for (int i = 0; i < nact; i++)
      {
        real h = 1.0 - v[i];
        m0 += h;
        mx += x[i][0] * h;
        my += x[i][1] * h;
        mz += x[i][2] * h;
      }

  PetscPrintf (PETSC_COMM_WORLD, "Print moments for tests: \n");
  PetscPrintf (PETSC_COMM_WORLD, "m0 = %lf\n", m0);
  PetscPrintf (PETSC_COMM_WORLD, "mx = %lf\n", mx / m0);
  PetscPrintf (PETSC_COMM_WORLD, "my = %lf\n", my / m0);
  PetscPrintf (PETSC_COMM_WORLD, "mz = %lf\n", mz / m0);

}
