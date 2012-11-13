/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-potential.h"

/*
  Context to mask the explicit form of PETSC stuff.  The corresponding
  typedef struct Context Context is in the header file.
*/
struct Context {
  DA da;                /* array descriptor */
  Vec v;                /* ref to a vector storing potential values */
  PetscScalar ***v_;    /* v_[k][j][i] points to the real data */
  int counter;          /* counter indicating how many are left */
  int local_size; /* length of  vector and array, better  save it than
                     use VecGetSize() from time to time */
  real h[3];
  real interval[2];
  int i0, j0, k0;
  int ni, nj, nk;
};

/* n = N * j + i, return i and j */
static void divmod (const int n, const int N, int* j, int* i)
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
  while  iterating over  it.  It  may though  VecDestroy()  the vector
  right after calling this function as we will increment the refcount.
*/
Context* bgy3d_pot_create (DA da, const ProblemData *PD, Vec v)
{
  /* memory for context pointer */
  Context *s = malloc (sizeof *s);

  s->da = da;
  PetscObjectReference ((PetscObject) s->da); /* DADestroy() it! */

  /* Do not copy the input  vector, save a reference instead, but also
     increment the refcount. We will have to VecDestroy() it too: */
  s->v = v;
  PetscObjectReference ((PetscObject) s->v);

  /* From  now  on  v_[k][j][i]  can   be  used  to  access  a  vector
     element: */
  DAVecGetArray (da, s->v, &s->v_);

  /* Get local portion of the grid */
  DAGetCorners (da, &s->i0, &s->j0, &s->k0, &s->ni, &s->nj, &s->nk);

  /* Set counter number and save redundant vector length: */
  s->local_size = s->ni * s->nj * s->nk;
  s->counter = s->local_size;

  /* copy N3 and h3 from PD */
  FOR_DIM
    s->h[dim] = PD->h[dim];

  /* copy interval */
  s->interval[0] = PD->interval[0];
  s->interval[1] = PD->interval[1];

  return s;
}

/* Value fetch interface */
int bgy3d_pot_get_value (Context *s, int n, real x[n][3], real v[n])
{
  /* how many elements would be skipped */
  const int nstart = s->local_size - s->counter;

  /* number of values that would be actually fetched */
  const int nact = MIN(n, s->counter);

  /* generating coordinates */
  int p = 0;
  for (int ijk = nstart; ijk < nstart + nact; ijk++)
    {
      int ij, i, j, k;

      /*
        Get coordinatens i, j, and k within the local block:

        ijk = k * NI * NJ + ij
        ij = j * NI + i
      */
      divmod (ijk, s->ni * s->nj, &k, &ij);
      divmod (ij, s->ni, &j, &i);

      /* Add corner coordinates: */
      i += s->i0;
      j += s->j0;
      k += s->k0;

      x[p][0] = i;
      x[p][1] = j;
      x[p][2] = k;
      v[p] = s->v_[k][j][i];
      // /* apply 'real' coordinates */
      // x[i][0] = i1 * s->h[0] + s->interval[0];
      // x[i][1] = j1 * s->h[1] + s->interval[0];
      // x[i][2] = k1 * s->h[2] + s->interval[0];
      p++;
    }
  assert (p == nact);

  /* update counter */
  s->counter -= nact;

  /* reset counter to original point once we fetched all the values */
  if (nact == 0)
    s->counter = s->local_size;

  return nact;
}

/* clean up the memory for public *pcontext */
void bgy3d_pot_destroy (Context *s)
{
  /* Required complement of DAVecGetArray(): */
  DAVecRestoreArray (s->da, s->v, &s->v_);

  /* Decrement  refcounts and  destroy  them if  the refcount  reached
     zero: */
  DADestroy (s->da);
  VecDestroy (s->v);

  /* free the whole context */
  free (s);
}

/*
 * Test for interface, to call this test function:
 *
 * 1. Create the context pointer:
 *
 *    Context *pcontext = bgy3d_pot_create (da, PD, vec);
 *
 * 2. Call the test function to check whether it prints out the right
 *    results:
 *
 *     bgy3d_pot_test (pcontext);
 *
 * 3. Remember to clean the memory:
 *
 *    bgy3d_pot_destroy (pcontext);
 */
void bgy3d_pot_test (Context *s)
{
  const int chunk_size = 120;
  real v[chunk_size];
  real x[chunk_size][3];
  int nact;

  /* Calculate moments for tests. First initializing local sums: */
  real m0 = 0.0;
  real m1[3] = {0.0, 0.0, 0.0};
  while ((nact = bgy3d_pot_get_value (s, chunk_size, x, v)))
    for (int i = 0; i < nact; i++)
      {
        real h = 1.0 - v[i];
        m0 += h;
        for (int j = 0; j < 3; j++)
          m1[j] += x[i][j] * h;
      }

  /* broadcast results of each worker to total sums */
  {
    int err;
    err = MPI_Allreduce (MPI_IN_PLACE, &m0, 1 , MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    assert (err == MPI_SUCCESS);

    err = MPI_Allreduce (MPI_IN_PLACE, m1, 3 , MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    assert (err == MPI_SUCCESS);
  }

  PetscPrintf (PETSC_COMM_WORLD, "Print moments for tests: \n");
  PetscPrintf (PETSC_COMM_WORLD, "m0 = %lf\n", m0);
  PetscPrintf (PETSC_COMM_WORLD, "mx = %lf\n", m1[0] / m0);
  PetscPrintf (PETSC_COMM_WORLD, "my = %lf\n", m1[1] / m0);
  PetscPrintf (PETSC_COMM_WORLD, "mz = %lf\n", m1[2] / m0);

}
