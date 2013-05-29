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
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-vec.h"          /* vec_ref(), ... */
#include "bgy3d-fftw.h"         /* bgy3d_fft_interp() */
#include "bgy3d-potential.h"
#include "bgy3d-poisson.h"      /* bgy3d_poisson() */

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
  real L[3];            /* box size */
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
  while  iterating over it.   It may  though vec_destroy()  the vector
  right after calling this function as we will increment the refcount.
*/
static Context* bgy3d_pot_create (const State *BHD, Vec v)
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
    Uincrement the refcount. We will have to vec_destroy() it
    too:
  */
  s->v = vec_ref (v);

  /* The first  time interpolation is requested we  put here something
     more usefull: */
  s->v_fft = NULL;

  /* From  now  on  v_[k][j][i]  can   be  used  to  access  a  vector
     element: */
  DMDAVecGetArray (s->da, s->v, &s->v_);

  /* Get local portion of the grid */
  DMDAGetCorners (s->da, &s->i0, &s->j0, &s->k0, &s->ni, &s->nj, &s->nk);

  /* copy N3 and h3 from PD */
  FOR_DIM
    s->h[dim] = BHD->PD->h[dim];

  /* copy interval */
  FOR_DIM
    s->L[dim] = BHD->PD->L[dim];

  /* Initalize counter: */
  s->ijk = 0;

  return s;
}

void bgy3d_pot_interp (Context *s, int n, /* const */ real x[n][3], real v[n])
{
  /* Prepare Fourier coefficients, if not has been already done: */
  if (s->v_fft == NULL)
    {
      s->v_fft = vec_create (s->dc); /* complex */
      MatMult (s->fft_mat, s->v, s->v_fft);
    }

  /* Translate site coordinates into real grid coordinates where
     integer values correspond to the grid nodes: */
  real y[n][3];
  for (int i = 0; i < n; i++)
    FOR_DIM
      y[i][dim] = (x[i][dim] + s->L[dim] / 2) / s->h[dim];

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

      x[p][0] = i * s->h[0] - s->L[0] / 2; /* x */
      x[p][1] = j * s->h[1] - s->L[1] / 2; /* y */
      x[p][2] = k * s->h[2] - s->L[2] / 2; /* z */
      v[p] = s->v_[k][j][i] * dV;          /* value * weight */

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
  /* Required complement of DMDAVecGetArray(): */
  DMDAVecRestoreArray (s->da, s->v, &s->v_);

  /* Decrement  refcounts and  destroy  them if  the refcount  reached
     zero: */
  DMDestroy (&s->da);
  DMDestroy (&s->dc);
  MatDestroy (&s->fft_mat);
  vec_destroy (&s->v);

  /* Only if interpolation was really used: */
  if (s->v_fft)
    vec_destroy (&s->v_fft);

  /* free the whole context */
  free (s);
}

/* Dipole of the cores of a single specie: */
static void dipole (int n, const Site sites[n], real d[3], real *d_norm)
{
  FOR_DIM
    d[dim] = 0.0;

  for (int i = 0; i < n; i++)
    FOR_DIM
      d[dim] += sites[i].charge * sites[i].x[dim];

  *d_norm = sqrt (SQR (d[0]) + SQR (d[1]) + SQR (d[2]));
}

static void moments (const State *BHD, Vec v,
                     real *q, real *x, real *y, real *z)
{
  const real *h = BHD->PD->h;   /* h[3] */
  const real h3 = h[0] * h[1] * h[2];

  /* 0th moment: */
  *q = h3 * vec_sum (v);

  /* 1st moments: */
  bgy3d_vec_moments1 (BHD->da, v, x, y, z);

  *x *= h3 * h[0];
  *y *= h3 * h[1];
  *z *= h3 * h[2];
}

/*
  Function to  generate solvent field that  is supposed to  act on the
  solute  electrons. At  the moment  only the  electrostatic  field is
  computed.   This  "reaction" field  should  be  consistent with  the
  "action" of  the solute electrons  on the solvent sites  computed in
  bgy3d-solvent.c  that  is also  pure  electrostatics. Returns  both,
  electrostatic potential  with a  boundary correction and  the actual
  solvent charge density.
*/
static void bgy3d_solvent_field (const State *BHD, /* intent(in) */
                                 int m, const Site solvent[m],
                                 Vec g[m],        /* intent(in) */
                                 Vec ve, Vec rho) /* intent(out) */
{
  /* FIXME: this  code assumes  the same density  rho for  all solvent
     particles. */

  /*
    Solvent   charge  density   for  Poisson   equation   and  (later)
    integration   with   solute   electrostatic  potential   will   be
    accumulated here:
  */
  VecSet (rho, 0.0);
  for (int i = 0; i < m; i++)
    VecAXPY (rho, solvent[i].charge * BHD->PD->rho, g[i]);

  /*
    Solve Poisson equation for rho.  Note that the output potential is
    in kcals  (see the  definiton of EPSILON0INV)  as all  energies in
    this code are:
  */
  bgy3d_poisson (BHD, ve, rho, -4 * M_PI * EPSILON0INV);

  /*
   Solving Poisson  equation by  FFT results in  a potential  with the
   mean value zero over the whole volume. That is how the ve(k) is set
   for k = 0. The following code adds a correction to that field which
   ensures  that  the  potential   is  zero  at  the  boundary.   This
   corresponds  to  a  boundary  as  a grounded  metallic  cage.   The
   corrected potential  is, in effect, a superposition  of the solvent
   electrostatic field and a surface charge on that metallic cage:
  */
  {
    local Vec x = vec_create (BHD->da);

    /*
      Correction is  a linear function  of potential on  the boundary.
      Almost  the  same  can  be  achieved by  invoking  the  function
      bgy3d_impose_laplace_boundary  (BHD,  ve,  x) except  that  that
      function also zeroes the boundary explicitly:
    */
    MatMult (BHD->dirichlet_mat, ve, x);

    /* Subtract it, making ve vanish at the boundary: */
    VecAXPY (ve, -1.0, x);

    vec_destroy (&x);
  }

  bgy3d_vec_save ("ve.bin", ve); /* for debugging only */
}


/* Print  a table  with some  site  info and  site-specific values  of
   potential from vs[n] in the last column: */
static void print_table (int n, const Site sites[n], const real vs[n])
{
  PetscPrintf (PETSC_COMM_WORLD,
               "#\t site\t x        \t y        \t z        \t q        \t δv\n");
  for (int i = 0; i < n; i ++)
    PetscPrintf (PETSC_COMM_WORLD,
                 "%d\t%5s\t% f\t% f\t% f\t% f\t% f\n",
                 i + 1, sites[i].name,
                 sites[i].x[0], sites[i].x[1],  sites[i].x[2],
                 sites[i].charge, vs[i]);
}


/*
  Prints  and returns  some info  interesting to  the caller.   At the
  moment  it  returns  the  iterator  over  the  medium  electrostatic
  potential. See bgy3d-potentialc:
*/
Context* info (const State *BHD,
               int m, const Site solvent[m],
               int n, const Site solute[n],
               Vec g[m],           /* in */
               Vec uc, Vec uc_rho) /* in */
{
  /* Solvent electrostatic potential field: */
  local Vec ve = vec_create (BHD->da);

  /* Keep solvent charge density for integration: */
  local Vec ve_rho = vec_create (BHD->da);

  /* This fills Vec ve with solvent electrostatic potential: */
  bgy3d_solvent_field (BHD, m, solvent, g, ve, ve_rho);

  /* Return the iterator over the solvent field: */
  Context *v = bgy3d_pot_create (BHD, ve);

  /* Also compute dipole moment: */
  real d[3], d_norm;
  dipole (n, solute, d, &d_norm);

  PetscPrintf (PETSC_COMM_WORLD,
               "|<x|ρ_n>| = |% f, % f, % f| = %f (dipole moment of solute cores)\n",
               d[0], d[1], d[2], d_norm);

  {
    real q, x, y, z;
    moments (BHD, ve_rho, &q, &x, &y, &z);

    const real d = sqrt (SQR (x) + SQR (y) + SQR (z));

    PetscPrintf (PETSC_COMM_WORLD,
                 "|<x|ρ_v>| = |% f, % f, % f| = %f (dipole moment of solvent medium)\n",
                 x, y, z, d);
    PetscPrintf (PETSC_COMM_WORLD, "q_v = %f (charge of solvent medium)\n", q);
  }
  /*
    Integration of:

    1.  Interaction  energy of solute  point cores with  the solvent
    electrostatic field.

    2. Solvent electrostatic field with diffuse solute electron
    density.

    3. Solvent charge density with long-range solute electrostatic
    field.

    Vec uc and uc_rho  are returned by bgy3d_solute_field().  Vec ve
    and ve_rho are obtained from bgy3d_solvent_field().
  */

  /* 1. */
  {
    real xs[n][3], vs[n];   /* coordinates and potential values */

    /* Extract coordinates into plain array: */
    for (int i = 0; i < n; i++)
      FOR_DIM
        xs[i][dim] = solute[i].x[dim];

    bgy3d_pot_interp (v, n, xs, vs);

    print_table (n, solute, vs);

    real val1 = 0.0;
    for (int i = 0; i < n; i++)
      val1 += solute[i].charge * vs[i];

    PetscPrintf (PETSC_COMM_WORLD,
                 "<U_v|ρ_N> = %lf (solvent electrostatic field with solute point nuclei)\n",
                 val1);
  }

  /* 2 and 3. */
  real val2, val3;
  {
    /* Dot-product is an integral over space (up to a factor): */
    const real h3 = BHD->PD->h[0] * BHD->PD->h[1] * BHD->PD->h[2];
    val2 = h3 * vec_dot (ve, uc_rho);
    val3 = h3 * vec_dot (uc, ve_rho);
  }

  PetscPrintf (PETSC_COMM_WORLD,
               "<U_v|ρ_u> = %lf "
               "(solvent electrostatic field with diffuse charge density of solute)\n",
               val2);
  PetscPrintf (PETSC_COMM_WORLD,
               "<ρ_v|U_u> = %lf "
               "(solvent charge density with long-range electrostatic field of solute)\n",
               val3);
  PetscPrintf (PETSC_COMM_WORLD,
               "     diff = % lf\n", val3 - val2);

  vec_destroy (&ve);            /* yes, we do! */
  vec_destroy (&ve_rho);

  return v;
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

      const real V = BHD->PD->L[0] * BHD->PD->L[1] * BHD->PD->L[2];
      PetscPrintf (PETSC_COMM_WORLD, "Moments divided by cell volume V = %lf: \n", V);
      PetscPrintf (PETSC_COMM_WORLD, "<1 * v> = %lf\n", m0 / V);
      PetscPrintf (PETSC_COMM_WORLD, "<x * v> = %lf\n", m1[0] / V);
      PetscPrintf (PETSC_COMM_WORLD, "<y * v> = %lf\n", m1[1] / V);
      PetscPrintf (PETSC_COMM_WORLD, "<z * v> = %lf\n", m1[2] / V);
    }

  bgy3d_pot_destroy (s);
}
