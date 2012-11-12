/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-potential.h"

/*
  Context to mask the explicit form of PETSC stuff.  The corresponding
  typedef struct Context Context is in the header file.
*/
struct Context {
  Vec v;                /* ref to a vector storing potential values */
  int counter;          /* counter indicating how many are left */
  int nmax; /* length  of vector and  array, better  save it  than use
               VecGetSize() from time to time */
  int N[3];
  real h[3];
  real interval[2];
};

static void rectangle (const int n, const int N, int* j, int* i);

/*
  Put the memory  allocation here.  This has be to  called from C side
  since  we don't  want to  allocate  memory for  vector from  fortran
  return the pointer to context after memory allocation.

  Vec v is saved in the returned iterator context by reference, not by
  copying.  User  code should  not change the  contents of  the vector
  while  iterating over  it.  It  may though  VecDestroy()  the vector
  right after calling this function as we will increment the refcount.
*/
Context* bgy3d_pot_create (DA da, const ProblemData *PD, Vec v)
{
  int i0, j0, k0;
  int ni, nj, nk;

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);
  const int m = ni * nj * nk;

  /* memory for context pointer */
  Context *pcontext = malloc(sizeof *pcontext);

  /* Do not copy the input  vector, save a reference instead, but also
     increment the refcount. We will have to VecDestroy() it too: */
  pcontext->v = v;
  PetscObjectReference ((PetscObject) pcontext->v);

  /* set counter number and save vector length */
  pcontext->counter = m;
  pcontext->nmax = m;

  /* copy N3 and h3 from PD */
  FOR_DIM
    {
      pcontext->N[dim] = PD->N[dim];
      pcontext->h[dim] = PD->h[dim];
    }

  /* copy interval */
  pcontext->interval[0] = PD->interval[0];
  pcontext->interval[1] = PD->interval[1];

  return pcontext;
}

/* Value fetch interface */
int bgy3d_pot_get_value (Context *s, int n, real x[n][3], real v[n])
{
  /* head index for context->v and context->x */
  const int nstart = s->nmax - s->counter;

  /* number of values that would be actually fetched */
  const int nact = MIN(n, s->counter);

  /* array storing indices that would be fetched from context->v */
  int idx[nact];

  /* generating coordinates */
  for (int i = 0; i < nact; i++)
    {
      int ijk = nstart + i;
      int ij, ix, iy, iz;

      /* ijk = iz * Nx * Ny + ij
       * ij = iy * Nx + ix */
      rectangle (ijk, s->N[0] * s->N[1], &iz, &ij);
      rectangle (ij, s->N[0], &iy, &ix);

      /* apply 'real' coordinates */
      x[i][0] = ix * s->h[0] + s->interval[0];
      x[i][1] = iy * s->h[1] + s->interval[0];
      x[i][2] = iz * s->h[2] + s->interval[0];

      /* save ijk to indicies array for usage later */
      idx[i] = ijk;
    }


  /* Get values from context->v[nstart : nstart + nact] */
  PetscErrorCode ierr = VecGetValues(s->v, nact, idx, v);
  assert (!ierr);

  /* update counter */
  s->counter -= nact;

  /* reset counter to original point once we fetched all the values */
  if (nact == 0)
    s->counter = s->nmax;

  return nact;
}

/* clean up the memory for public *pcontext */
void bgy3d_pot_destroy (Context *s)
{
  /* Decrement refcount and destroy the vector if the refcount reached
     zero: */
  VecDestroy (s->v);

  /* free the whole context */
  free (s);
}

/* n = N * j + i, return i and j */
static void rectangle (const int n, const int N, int* j, int* i)
{
  *j = n / N;
  *i = n % N;
}

/* Test for interface */
void bgy3d_pot_test (Context *s)
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
