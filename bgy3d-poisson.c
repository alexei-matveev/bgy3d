/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O_solutes.c,v 1.3 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

/*
  There are  as many ways to  solve the Poisson equation  as there are
  ways  to construct  a discrete  representation of  Laplace operator.
  The FFTW  convention for the  inverse FFT expressing the  real space
  u(x) via the momentum space ũ(k) is:

                kx              2πi/n
    u  = Σ  ũ  ω  ,   with ω = e
     x    k  k

  where the summation is over 0 <= k < n. Thus, the forward difference

                      kx   k
    u   - u  = Σ  ũ  ω   (ω - 1)
     x+1   x    k  k

  and the backward difference is

                      kx       -k
    u - u    = Σ  ũ  ω   (1 - ω  ).
     x   x-1    k  k

  The second-order difference is thus

                              kx   k        -k
    u   - 2u  + u    = Σ  ũ  ω   (ω  - 2 + ω  )
     x+1    x    x-1    k  k

                                  kx    2
                     = - 4 Σ  ũ  ω   sin (πk/n).
                            k  k

  Note that a plane wave  ω^kx is an eigenfunction of the second-order
  difference operator with an eigenvalue

    2 cos(2πk/n) - 2 =  -4 sin^2 (πk/n).

  The eigenvalues for k and k' = n - k are equal. For small ratios k/n
  the eigenvalues  are indeed  approximately proportional to  k^2. For
  the higher-order O(h^4) stencil

    (- u    + 16u    - 30u  + 16u    - u   ) / 12
        x+2      x+1      x      x-1    x-2

  The corresponding eigenvalues would have been

    (-2 cos(4πk/n) + 32 cos(2πk/n) - 30) / 12

  which has  a similar quadratic shape  for small k/n but  a lower and
  sharper minimum at around k  = n/2.  This similarity is probably the
  reason the original code uses the spectral representation of Laplace
  operator with eigenvalues proportional to k^2 for k <= n/2 and to (n
  - k)^2 otherwise thus  having a "casp" at k =  n/2.  Such a spectrum
  will correspond to a non-compact stencil in real space. Here are the
  first 6 elements of this (even) stencil:

    -3.289868 2.000000 -0.500000 0.222222 -0.125000 0.080000 ...

  Computed by the octave code for a large even n:

    a = [0:n/2, n/2-1:-1:1];
    b = real (ifft (-4 * ((pi / n) * a) .^ 2));

  What are the arguments in favor of and against a particular approach
  is not quite clear.
*/

#include <stdbool.h>
#include "bgy3d.h"
#include "bgy3d-poisson.h"
#include <complex.h>            /* after fftw.h */

/* GCC extensions: */
#if __GNUC__ >= 3
#define likely(x)    __builtin_expect (!!(x), 1)
#define unlikely(x)  __builtin_expect (!!(x), 0)
#else
#define likely(x)    (x)
#define unlikely(x)  (x)
#endif

/* The alternative  works as well,  though the regression test  do not
   survive with more thatn 3-4 digits: */
#ifndef POISSON_AS_INVERSE_LAPLACE
/*
  Solve  Poisson  Equation  in  Fourier space  and  get  elestrostatic
  potential by inverse FFT.

  Vec uc is intent(out).
  Vec rho is intent(in).
  real q is the overall factor.

  To get the potential in kcal/mol as used in the rest of the code you
  need to supply q = -4π/ε₀ that is -4 * M_PI * EPSILON0INV.

  As a  matter of  fact, it  appears that one  could provide  the same
  factual parameter  for rho and  uc to effectively solve  the Poisson
  equation "in place".

  Except of  temporary allocation of a  complex Vec does  not have any
  side effect.
*/
void bgy3d_poisson (const State *BHD, Vec uc, Vec rho, real q)
{
  const real *interval = BHD->PD->interval; /* [2] */
  const int *N = BHD->PD->N;    /* N[3] */

  /* Otherwise needs some work: */
  assert (N[0] == N[1]);
  assert (N[0] == N[2]);

  const real L = interval[1] - interval[0];
  const int NNN = N[0] * N[1] * N[2];

  /* Scratch complex vector: */
  Vec work;
  DACreateGlobalVector (BHD->dc, &work);

  /* Get FFT of  rho: rho(i, j, k) -> fft_rho(kx,  ky, kz) placed into
     complex work: */
  MatMult (BHD->fft_mat, rho, work);

  /*
    Solving Poisson Equation (note the  absence of -4π factor) with FFT
    and IFFT:

        Δu(x, y, z) = ρ(x, y, z)

    because of x = ih, y = jh, and z = kh, with grid spacing h = L/n:

        n² / L²  Δu(i, j, k) = ρ(i, j, k)

    In Fourier  space the relation  between FFT images  of ρ and  u is
    (see FFTW manual "What FFTW Really Computes"):

        u(kx, ky, kz) = ρ(kx, ky, kz) / (4 π² k² / L²)

    with

        k² = kx² + ky² + kz²

    being the  sum of  squared integers.  Finally  do the  inverse FFT
    (see  FFTW manual "What  FFTW Really  Computes").  Because  of the
    normalization IFFT(FFT(f)) = n³ * f we have:

        u(i, j, k) = 1 / n³  * IFFT(u(kx, ky, kz))
  */

  /* With q = -4π/ε₀ you would get the potential: */
  const real scale = - q / NNN / 4;

  /* Loop over local portion of the k-grid */
  {
    int x[3], n[3], i[3], ic[3];
    DAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

    complex ***work_;
    DAVecGetArray (BHD->dc, work, &work_);

    for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
      for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
        for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
          {
            /* FIXME: what if we  change the complex vectors to remove
               the redundnacy? */
            FOR_DIM
              {
                if (i[dim] <= N[dim] / 2)
                  ic[dim] = i[dim];
                else
                  ic[dim] = i[dim] - N[dim];
              }

            if (ic[0] == 0 && ic[1] == 0 && ic[2] == 0)
              {
                /* The gamma point, k = 0, we cannot divide by 0: */
                work_[i[2]][i[1]][i[0]] = 0.0; /* complex */
              }
            else
              {
                /* For  i, j,  and k  less than  or equal  to  N/2 and
                   uniform box of size  L this expression evaluates to
                   (π/L)² (i² + j² + k²) */
                const real k2 = SQR (M_PI / L) *
                  (SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]));

                const real fac = scale / k2;

                /* Here we compute in place: uc(kx, ky, kz) := scale *
                   rho(kx, ky, kz) / k^2 */
                work_[i[2]][i[1]][i[0]] *= fac; /* complex */
              }
          }
    DAVecRestoreArray (BHD->dc, work, &work_);
  }

  /* u(x, y, z) := IFFT(u(kx, ky, kz)) */
  MatMultTranspose (BHD->fft_mat, work, uc);

  VecDestroy (work);
}
#else
/*
  This solves the  equation Δu = qρ with  Laplacian defined by 7-point
  stencil and  periodic boundary  conditions.  Thus with  q = 1  it is
  exactly the pseudo-inverse of  the Laplacian operator implemented by
  7-point  stencil.   Given  the  periodicity the  pseudo-inverse  may
  differ by at most a constant.

  To get the potential in kcal/mol as used in the rest of the code you
  need to supply q = -4π/ε₀ that is -4 * M_PI * EPSILON0INV.
*/
void bgy3d_poisson (const State *BHD, Vec uc, Vec rho, real q)
{
  const real *h = BHD->PD->h;   /* h[3] */
  const int *N = BHD->PD->N;    /* N[3] */

  const int NNN = N[0] * N[1] * N[2];

  /* Scratch complex vector: */
  Vec work;
  DACreateGlobalVector (BHD->dc, &work);

  /* Get FFT of  rho: rho(i, j, k) -> fft_rho(kx,  ky, kz) placed into
     complex work: */
  MatMult (BHD->fft_mat, rho, work);

  /* With q = -4π/ε₀ you would get the potential: */
  const real scale = - q / NNN / 4;

  /* Loop over local portion of the k-grid */
  {
    int i0, j0, k0, ni, nj, nk;
    DAGetCorners (BHD->dc, &i0, &j0, &k0, &ni, &nj, &nk);

    complex ***work_;
    DAVecGetArray (BHD->dc, work, &work_);

    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          {
            /*
              For small  i, j, and k  and uniform box of  size L this
              expression approximates to (π/L)² (i² + j² + k²).

              FIXME: SQR() macro evaluates the argument twice!
            */
            const real k2 =                             \
              SQR (sin (M_PI * i / N[0]) / h[0]) +
              SQR (sin (M_PI * j / N[1]) / h[1]) +
              SQR (sin (M_PI * k / N[2]) / h[2]);

            real fac;
            if (likely (k2 != 0.0))
              fac = scale / k2;
            else
              fac = 0.0;        /* gamma-point */

            /* Here we compute in place: u(kx, ky, kz) := scale *
               rho(kx, ky, kz) / k^2 */
            work_[k][j][i] *= fac; /* complex */
          }
    DAVecRestoreArray (BHD->dc, work, &work_);
  }

  /* u(x, y, z) := IFFT(u(kx, ky, kz)) */
  MatMultTranspose (BHD->fft_mat, work, uc);

  VecDestroy (work);
}
#endif

#ifdef L_BOUNDARY

typedef struct Boundary {
  int border;
  const int *N;                 /* N[3] => PD->N[3] */
} Boundary;

/* Returns true iff the point (i, j, k) is inside the volume described
   by the boundary: */
static bool inside_boundary (const Boundary *b, int i, int j, int k)
{
  /*
    With current boundary definition return true iff:

    border < i < N[0] - border
    border < j < N[1] - border
    border < k < N[2] - border
  */
  const int *N = b->N;          /* N[3] */
  const int border = b->border;

  return \
    likely (border < i) && likely (i < N[0] - border) &&
    likely (border < j) && likely (j < N[1] - border) &&
    likely (border < k) && likely (k < N[2] - border);
}

/* Returns   a   description   of    the   boundary   for   use   with
   inside_boundary(): */
static Boundary make_boundary (const ProblemData *PD)
{
  const int *N = PD->N;         /* N[3] */
  const real *h = PD->h;        /* h[3] */
  const real L = PD->interval[1] - PD->interval[0];
  const real zpad = PD->zpad;

  /*
    FIXME:  With  this definition  and  the  default  zpad =  L/2  the
    variable "border" comes out equal 1.  For the actual boundary, see
    inside_boundary(), this means the planes with i = 0, i = 1 and i =
    N - 1 are  all part of the boundary (the same  holds for for j and
    k,  of  course). Given  the  periodicity,  one  can think  of  the
    boundary as the planes i = -1, 0, 1.  Thus the default boundary is
    three layers "thick"! I am not sure that this was the intention.

    Now that the boundary definition is localized one could switch the
    default boundary to just one layer by omitting +1 in the following
    expression. The code appears  to still work, though the regression
    tests fail though with deviations  in the third or fourth digit of
    the moments.
  */
  const int border = 1 + (int) ceil ((L - 2.0 * zpad) / h[0] / 2.0);

  /* Holds for all regression tests! */
  assert (border == 1);

  const Boundary boundary = {border, N};
  return boundary;
}


/* Copies the values of Vec g at a boundary to Vec b.  The rest of Vec
   b is NOT changed! */
static void copy_boundary (DA da, const Boundary *vol, Vec g, Vec b)
{
  /* Get local portion of the grid */
  int i0, j0, k0, ni, nj, nk;
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  PetscScalar ***g_, ***b_;
  DAVecGetArray (da, g, &g_);
  DAVecGetArray (da, b, &b_);

  /* loop over local portion of grid */
  for (int k = k0; k < k0 + nk; k++)
    for (int j = j0; j < j0 + nj; j++)
      for (int i = i0; i < i0 + ni; i++)
        if (unlikely (!inside_boundary (vol, i, j, k)))
          b_[k][j][i] = g_[k][j][i];

  DAVecRestoreArray (da, g, &g_);
  DAVecRestoreArray (da, b, &b_);
}

/*
  This function sets everything in  the 3d-array g[:, :, :] outside of
  the central section g[i, j, k] with b  < i < N - b, (same for j, and
  k) to  the value "value" (typically  0.0, less often  1.0). With the
  value of zpad equal  to the (half) the box size L  the value of b is
  1.
*/
static void set_boundary (DA da, const Boundary *vol, Vec g, real value)
{
  /* Get local portion of the grid */
  int i0, j0, k0, ni, nj, nk;
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  PetscScalar ***g_;
  DAVecGetArray (da, g, &g_);

  /* loop over local portion of grid */
  for (int k = k0; k < k0 + nk; k++)
    for (int j = j0; j < j0 + nj; j++)
      for (int i = i0; i < i0 + ni; i++)
        if (unlikely (!inside_boundary (vol, i, j, k)))
          g_[k][j][i] = value;

  DAVecRestoreArray (da, g, &g_);
}

/*
  The  alternative works  on a  single CPU  too. Again  the regression
  tests  "fail"  with 3-4  digits  the  same  as before.   The  bigger
  difference  appears  to be  that  the matrix  in  AIJ  format has  a
  dedicated  solver for  linear equations.   The  algorithm KSPSolve()
  uses with matrix free implementation  appears to be slower and fails
  for more than one CPU.
*/
#ifndef MATRIX_FREE_LAPLACE
/*
  Create    and   initialize    Laplace   matrix    with   appropriate
  stencil.  FIXME: this  should be  probably implemented  using matrix
  free facilities of Petsc.
*/
static void lap_mat_create (const DA da, const real h[3],
                            const Boundary *const vol,
                            Mat *M) /* intent(out) */
{
  const PetscScalar one = 1.0;

  PetscPrintf (PETSC_COMM_WORLD, "Assembling Matrix...");

  /* Create Matrix with appropriate non-zero structure */
  DAGetMatrix (da, MATMPIAIJ, M);

  MatZeroEntries (*M);

  /*
    This code constructs  (a compact representation of) the  N^3 x N^3
    matrix M.  An identity or another diagonal matrix will always have
    row == col.  The Laplace matrix is sparse as in "nearly diagonal".
    col[3]  will be  used to  hold  the coordinated  of three  stencil
    points in  a row.  Note that  "3" is not  the space  dimension but
    rather the  stencil size for  second derivative in  one dimension.
    The central  point, in  col[1] == *(col  + 1), corresponds  to the
    current grid point with indices i[].
  */

  /* Get local portion of the grid */
  int x[3], n[3], N[3];
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* Get grid dimensions N[3], sanity checks: */
  {
    int dim, dof, sw;
    DAPeriodicType wrap;
    DAStencilType st;
    DAGetInfo (da, &dim,
               &N[0], &N[1], &N[2], /* need this, rest for checks */
               NULL, NULL, NULL,
               &dof, &sw, &wrap, &st);

    /* It may  or may  not work  for other settings  too, it  was only
       tested with these: */
    assert (dim == 3);               /* 3d Vec */
    assert (dof == 1);               /* degrees of freedom */
    assert (sw >= 1);                /* stencil width */
    assert (st == DA_STENCIL_STAR);  /* stencil type */
    assert (wrap == DA_XYZPERIODIC); /* periodicity */
  }

  /* Loop over local portion of grid: */
  for (int k = x[2]; k < x[2] + n[2]; k++)
    for (int j = x[1]; j < x[1] + n[1]; j++)
      for (int i = x[0]; i < x[0] + n[0]; i++)
        {
          MatStencil row;
          /* This  is the  point  on the  diagonal  of the  N^3 x  N^3
             matrix: */
          row.i = i;
          row.j = j;
          row.k = k;

          /* Boundary,  if not  NULL. Otherwise  it will  be  a sparse
             representation of the plain Laplace operator: */
          if (vol && unlikely (!inside_boundary (vol, i, j, k)))
            {
              /* This  sets this  particular diagonal  element  of the
                 matrix to 1.0: */
              MatSetValuesStencil (*M, 1, &row, 1, &row, &one, ADD_VALUES);
            }
          else
            {
              MatStencil col[3]; /* Three stencil points in the row. */

              /* Central point of the stencil, same as row: */
              col[1].i = row.i;
              col[1].j = row.j;
              col[1].k = row.k;

              /*
                Other two  points of the  stencil offset by 1  in -dim
                and +dim, respectively are  set in the switch statment
                inside the loop over dim. The manual says:

                  The columns  and rows in  the stencil passed  in [to
                  MatSetValuesStencil()] MUST  be contained within the
                  ghost  region  of  the  given process  as  set  with
                  DMDACreateXXX() or MatSetStencil().

                See  asserts  above an  the  call  to DACreate3d()  in
                bgy3d-fftw.c,  especially the periodicity  and stencil
                options there.

                If you wrap the  indices modulo N, Petsc will complain
                that  they are not  in the  (local) range  in parallel
                runs with  NULL boundary.  With  NULL boundary indices
                may appear  out of range  and even negative.   I guess
                this is the intention.
              */
              FOR_DIM
                {
                  switch (dim)
                    {
                    case 0:
                      col[0].i = row.i - 1;
                      col[0].j = row.j;
                      col[0].k = row.k;

                      col[2].i = row.i + 1;
                      col[2].j = row.j;
                      col[2].k = row.k;
                      break;

                    case 1:
                      col[0].i = row.i;
                      col[0].j = row.j - 1;
                      col[0].k = row.k;

                      col[2].i = row.i;
                      col[2].j = row.j + 1;
                      col[2].k = row.k;
                      break;

                    case 2:
                      col[0].i = row.i;
                      col[0].j = row.j;
                      col[0].k = row.k - 1;

                      col[2].i = row.i;
                      col[2].j = row.j;
                      col[2].k = row.k + 1;
                      break;
                    }

                  /* Values to enter for the Laplacian stencil: */
                  const real h2 = SQR (h[dim]);
                  const PetscScalar v[3] = {1.0 / h2, -2.0 / h2, 1.0 / h2};

                  MatSetValuesStencil (*M, 1, &row, 3, col, v, ADD_VALUES);
                }
            }
        }

  MatAssemblyBegin (*M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd (*M, MAT_FINAL_ASSEMBLY);

  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}
#else
/*
  For a  plain laplacian only the  grid dimensions and  mesh sizes are
  necessary. For a Laplacian with boundary also the third field is set
  and used:
*/
typedef struct Lap {
  DA da;                        /* array descriptor */
  real h[3];                    /* grid spacing */
  Boundary vol;                 /* for boundary problem */
} Lap;


static Lap* context (Mat A)
{
  Lap *lap;
  MatShellGetContext (A, (void**) &lap);
  return lap;
}


/* Does y = A * x = Δx. Laplace operator by finite differences: */
static PetscErrorCode mat_mult_lap (Mat A, Vec x, Vec y)
{
  /* Only  matrices  constructed   by  mat_create_lap()  are  accepted
     here: */
  const Lap *lap = context (A);

  const real *h = lap->h;       /* h[3] */

  assert (h[0] == h[1]);
  assert (h[0] == h[2]);

  const real h2 = SQR (h[0]);

  /* loop over local portion of grid */
  {
    /*
      This will be a local  vector with ghost points.  Apparently the
      distinction between local and  global are the ghost points, and
      not what you may have thought.
    */
    Vec w;

    /*
      There is  a Pets infrastracture to get  temporary local vectors.
      Still  this may  lead to  allocations that  are only  freed upon
      destruction of the array descriptor:
    */
    DAGetLocalVector (lap->da, &w);

    /*
      This initiates  communication and  waits for its  completion ---
      ghost points  need to  be scattared to  the "neighbors".  Do not
      attempt to VecCopy (x, w), you naive Petsc user!
    */
    DAGlobalToLocalBegin (lap->da, x, INSERT_VALUES, w);
    DAGlobalToLocalEnd (lap->da, x, INSERT_VALUES, w);

    real ***w_, ***y_;
    DAVecGetArray (lap->da, w, &w_);
    DAVecGetArray (lap->da, y, &y_);

    /* Get local portion of the grid */
    int i0, j0, k0, ni, nj, nk;
    DAGetCorners (lap->da, &i0, &j0, &k0, &ni, &nj, &nk);

    // printf ("LOCAL: x = %d %d %d, n = %d %d %d\n", i0, j0, k0, ni, nj, nk);
    // DAGetGhostCorners (lap->da, &i0, &j0, &k0, &ni, &nj, &nk);
    // printf ("GHOST: x = %d %d %d, n = %d %d %d\n", i0, j0, k0, ni, nj, nk);

    /* This  may only  work with  periodic Petsc  vectors  and stencil
       width >= 1: */
    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          {
            /* One of 7 stencil points, the center: */
            const real w_kji = w_[k][j][i];

            /*
              The  indices into the  ghosted array  may appear  out of
              bounds   and   even    negative,   but   this   is   the
              intention. Parens,  though redundant mathematically, are
              intended to reduce the loss of precision:
            */
            y_[k][j][i] = (1.0 / h2) *          \
              (((w_[k][j][i + 1] - w_kji) +
                (w_[k][j][i - 1] - w_kji)) +
               ((w_[k][j + 1][i] - w_kji) +
                (w_[k][j - 1][i] - w_kji)) +
               ((w_[k + 1][j][i] - w_kji) +
                (w_[k - 1][j][i] - w_kji)));
          }
    DAVecRestoreArray (lap->da, w, &w_);
    DAVecRestoreArray (lap->da, y, &y_);

    /* The  counterpart to  DAGetLocalVector(), do  not  destroy, just
       give it back: */
    DARestoreLocalVector (lap->da, &w);
  }

  return 0;
}


/* This is what is actually used: */
static PetscErrorCode mat_mult_bnd (Mat A, Vec x, Vec y)
{
  mat_mult_lap (A, x, y);       /* y := Δx */

  Lap *lap = context (A);

  /* Untill  now we compute  plain-old Δx,  now the  boundary specific
     staff. Copy the values at the boundary from x to y as is: */
  copy_boundary (lap->da, &lap->vol, x, y);

  return 0;
}

static PetscErrorCode mat_destroy (Mat A)
{
  /* Only  matrices  constructed   by  lap_mat_create()  are  accepted
     here: */
  Lap *lap = context (A);

  DADestroy (lap->da);          /* destroys or decrements refcount */

  free (lap);

  return 0;
}

/* Creates a matix  shell, but does not associate  any operations with
   it: */
static void mat_create_shell (const DA da, void *ctx, Mat *A)
{
  /* Get dimensions and other vector properties: */
  int x[3], n[3], N[3];
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* Get grid dimensions N[3], sanity checks: */
  {
    int dim, dof, sw;
    DAPeriodicType wrap;
    DAStencilType st;
    DAGetInfo (da, &dim,
               &N[0], &N[1], &N[2], /* need this, rest for checks */
               NULL, NULL, NULL,
               &dof, &sw, &wrap, &st);

    /* It may  or may  not work  for other settings  too, it  was only
       tested with these: */
    assert (dim == 3);               /* 3d Vec */
    assert (dof == 1);               /* degrees of freedom */
    assert (sw >= 1);                /* stencil width */
    assert (st == DA_STENCIL_STAR);  /* stencil type */
    assert (wrap == DA_XYZPERIODIC); /* periodicity */
  }

  /*
    Create the matrix  itself. First, set total and  local size of the
    matrix (section).
  */
  const int N3 = N[2] * N[1] * N[0];
  const int n3 = n[2] * n[1] * n[0];

  /* Also puts internals (Lap struct) into the matrix: */
  MatCreateShell (PETSC_COMM_WORLD, n3, n3, N3, N3, ctx, A);
}


/* Creates Laplacian matrix: */
static void lap_mat_create (const DA da, const real h[3],
                            const Boundary *vol,
                            Mat *A)    /* intent(out) */
{
  /* Allocates storage for  a Lap struct. Make sure to  it is freed in
     mat_destroy(): */
  Lap *lap = malloc (sizeof *lap);

  /*
    That is all we need  for the Laplace operator, dimensions and mesh
    sizes.   Copy contents,  the  there are  no  guarantees about  the
    life-time of the object pointed to:
  */
  lap->da = da;
  for (int i = 0; i < 3; i ++)
    lap->h[i] = h[i];

  /* We will  DADestroy() it in  mat_destroy(), but it came  from user
     space, so increment the refcount: */
  PetscObjectReference ((PetscObject) lap->da);

  /* This is used to copy  over the boundary values: */
  if (vol)
    lap->vol = *vol;

  /* Create  matrix shell  with  proper dimensions  and associate  the
     context with it: */
  mat_create_shell (da, lap, A);

  /* Set matrix operations. If vol  is NULL, the resulting operator is
     plain Laplace: */
  if (vol)
    MatShellSetOperation (*A, MATOP_MULT,
                          (void (*)(void)) mat_mult_bnd);
  else
    MatShellSetOperation (*A, MATOP_MULT,
                          (void (*)(void)) mat_mult_lap);

  MatShellSetOperation (*A, MATOP_DESTROY,
                        (void (*)(void)) mat_destroy);
}
#endif

static void InitializeKSPSolver (Mat M, KSP *ksp)
{
  PC pc;

  /* Create ksp environment */
  KSPCreate (PETSC_COMM_WORLD, ksp);
  KSPGetPC (*ksp, &pc);

  /* FIXME: literal tolerances here: */
  KSPSetTolerances (*ksp, 1.0e-4, 1.0e-4, 1.0e+5, 1000);

  /* Set Matrix */
  //KSPSetOperators (*ksp, M, M, SAME_NONZERO_PATTERN);
  KSPSetOperators (*ksp, M, M, SAME_PRECONDITIONER);

  /* Set preconditioner */
  PCSetType (pc, PCBJACOBI);

  KSPSetInitialGuessNonzero (*ksp, PETSC_TRUE);

  /* runtime options will override default parameters */
  //KSPSetFromOptions(BHD->ksp);
}


/* Assemble Laplacian matrix and create KSP environment: */
void bgy3d_laplace_create (const DA da, const ProblemData *PD, Mat *M, KSP *ksp)
{
  const Boundary vol = make_boundary (PD);

  lap_mat_create (da, PD->h, &vol, M);

  InitializeKSPSolver (*M, ksp);
}

/*
  This solves the Dirichlet problem Δx  = 0, x(∂Ω) = v(∂Ω) in order to
  build v' =  v - x that  vanishes at the boundary: v'(∂Ω)  = 0.  More
  specifically, a linear equation for x is solved first:

    KSP * x = P * v

  an then v is being updated:
                         -1
    v := v - x = [1 - KSP   *  P] * v

  with KSP being  the Laplace boundary problem matrix  and P being the
  projector onto the boundary. That is b = P * v constructs a vector b
  that is equal to v at the boundary and zero everywhere else.

  Laplace equation is solved iteratively, so that, I assume, the input
  value for  x, does  matter. At the  later stages of  iterations when
  potential v is  only slightly changing at the  boundary preserving x
  across invocations should save some time. So is the theory.

  See tratment of the boundary  conditions in the thesis: pp.  116-177
  Eqs.  (5.107) - (5.110).

  Vec v, x are intent(inout), Vec b is a work array.
 */
void bgy3d_impose_laplace_boundary (const State *BHD, Vec v, Vec b, Vec x)
{
  const Boundary vol = make_boundary (BHD->PD);

  /*
    Get boundary b of v, the rest  of b is set to zero. Together it is
    a linear, albeit not invertible projection operation:

    b := P * v
  */
  VecSet (b, 0.0);
  copy_boundary (BHD->da, &vol, v, b);

  /*
    Solve Laplace  equation, update  x iteratively.  Ideally  from the
    state of the previous iteration:
            -1
    x := KSP   *  b
  */
  KSPSolve (BHD->ksp, b, x);

  /*
    If you think you need to  know how many iteration were required to
    solve the equation do this:

    int iter;
    KSPGetIterationNumber (BHD->ksp, &iter);
  */

  /*
    Subtract  solution  from  v.  The  whole is  in  effect  a  linear
    operation:
                  -1
    v := [1 - (KSP  *  P)] * v
  */
  VecAXPY (v, -1.0, x);

  /*
    FIXME:  formally  redundant.   Set  v  to  zero  at  the  boundary
    again.  Historically, the call  to bgy3d_impose_laplace_boundary()
    was  followed immediately  by an  equivalent of  set_boundary() at
    many places in the code.  So do this here instead. Commenting this
    call changes  results only slightly,  most probably due  to finite
    convergence thresholds in KSPSolve():
  */
  set_boundary (BHD->da, &vol, v, 0.0);

  /* WRONG? If  you preserve the  value of x  until the next  call the
     iterative solver will re-use  it as initial approximation for the
     next x. */
}
#endif


/* Laplace matrix, no boundary: */
void bgy3d_lap_mat_create (const DA da, const real h[3], Mat *A)
{
  lap_mat_create (da, h, NULL, A);
}
