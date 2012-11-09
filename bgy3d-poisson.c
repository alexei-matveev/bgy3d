/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O_solutes.c,v 1.3 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

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

/*
  Solve  Poisson  Equation  in  Fourier space  and  get  elestrostatic
  potential by inverse FFT.

  Vec uc is intent(out).
  Vec rho is intent(in).
  real q is the overall factor.

  As a  matter of  fact, it  appears that one  could provide  the same
  factual parameter  for rho and  uc to effectively solve  the Poisson
  equation "in place".

  Except of  temporary allocation of a  complex Vec does  not have any
  side effect.
*/
void bgy3d_poisson (const State *BHD, Vec uc, Vec rho, real q)
{
  const real *interval = BHD->PD->interval; /* [2] */
  const real *h = BHD->PD->h;   /* h[3] */
  const int *N = BHD->PD->N;    /* N[3] */

  const real L = interval[1] - interval[0];
  const real h3 = h[0] * h[1] * h[2];

  /* Scratch complex vector: */
  Vec work;
  DACreateGlobalVector (BHD->dc, &work);

  /* Get FFT of  rho: rho(i, j, k) -> fft_rho(kx,  ky, kz) placed into
     complex work: */
  MatMult (BHD->fft_mat, rho, work);

  /*
    Solving Poisson Equation (SI units) with FFT and IFFT:

        - Δu(x, y, z) = (1 / ε₀) ρ(x, y, z)

    because of x = ih, y = jh, and z = kh, with grid spacing h = L/n:

        - n² / L²  Δu(i, j, k) = (1 / ε₀) ρ(i, j, k)

    In Fourier  space the relation  between FFT images  of ρ and  u is
    (see FFTW manual "What FFTW Really Computes"):

        u(kx, ky, kz) = 1 / (4 π² ε₀ k² / L²) ρ(kx, ky, kz)

    with

        k² = kx² + ky² + kz²

    IFFT (see FFTW manual "What FFTW Really Computes"):

    because: IFFT(u(kx, ky, kz)) = n³ * u(i, j, k)

        u(i, j, k) = h³ / L³  * IFFT(u(kx, ky, kz))
  */

  /* EPSILON0INV = 1 / 4 π ε₀: */
  const real scale = q * EPSILON0INV / M_PI * h3 / (L * L * L);

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
                const real k2 = (SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0])) / SQR(L);

                const real fac = scale / k2;

                /* Here we compute in place: uc(kx, ky, kz) := scale *
                   rho(kx, ky, kz) / k^2 */
                work_[i[2]][i[1]][i[0]] *= fac; /* complex */
              }
          }
    DAVecRestoreArray (BHD->dc, work, &work_);
  }

  /* uc := IFFT(uc(kx, ky, kz)) */
  MatMultTranspose (BHD->fft_mat, work, uc);

  VecDestroy (work);
}


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


/* Returns a  non-negative number,  e.g. mod(-1, 10)  -> 9.   Does not
   work for b <= 0: */
static int mod (int a, int b)
{
  return ((a % b) + b) % b;
}


/*
  Create    and   initialize    Laplace   matrix    with   appropriate
  stencil.  FIXME: this  should be  probably implemented  using matrix
  free facilities of Petsc.
*/
static void InitializeLaplaceMatrix (const DA da, const real h[3],
                                     const Boundary *vol,
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

          /* Boundary */
          if (unlikely (!inside_boundary (vol, i, j, k)))
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
                inside  the loop  over dim.   FIXME: we  might  have a
                problem for parallel runs here as the manual says:

                  The columns  and rows in  the stencil passed  in [to
                  MatSetValuesStencil()] MUST  be contained within the
                  ghost  region  of  the  given process  as  set  with
                  DMDACreateXXX() or MatSetStencil().

                See  call to  DACreate3d() in  bgy3d-fft.c, especially
                the periodicity and stencil options there.
              */
              FOR_DIM
                {
                  switch (dim)
                    {
                    case 0:
                      col[0].i = mod (row.i - 1, N[0]);
                      col[0].j = row.j;
                      col[0].k = row.k;

                      col[2].i = mod (row.i + 1, N[0]);
                      col[2].j = row.j;
                      col[2].k = row.k;
                      break;

                    case 1:
                      col[0].i = row.i;
                      col[0].j = mod (row.j - 1, N[1]);
                      col[0].k = row.k;

                      col[2].i = row.i;
                      col[2].j = mod (row.j + 1, N[1]);
                      col[2].k = row.k;
                      break;

                    case 2:
                      col[0].i = row.i;
                      col[0].j = row.j;
                      col[0].k = mod (row.k - 1, N[2]);

                      col[2].i = row.i;
                      col[2].j = row.j;
                      col[2].k = mod (row.k + 1, N[2]);
                      break;
                    }

                  /* Sanity check: */
                  for (int p = 0; p < 3; p++)
                    {
                      assert (col[p].i >= 0);
                      assert (col[p].j >= 0);
                      assert (col[p].k >= 0);
                      assert (col[p].i < N[0]);
                      assert (col[p].j < N[1]);
                      assert (col[p].k < N[2]);
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

  InitializeLaplaceMatrix (da, PD->h, &vol, M);

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
    again. State BHD  is not modified by this  call. Historically, the
    call  to bgy3d_impose_laplace_boundary() was  followed immediately
    by bgy3d_boundary_set()  at many places  in the code.  So  do this
    here instead. Commenting this  call changes results only slightly,
    most probably due to finite convergence thresholds in KSPSolve():
  */
  bgy3d_boundary_set (BHD, v, 0.0);

  /* If you preserve the value of  x until the next call the iterative
     solver will re-use it as initial approximation for the next x. */
}
#endif

/*
  This function appears  to set everything in the  3d-array g[:, :, :]
  outside of the central section g[i, j, k]  with b < i < N - b, (same
  for  j, and  k)  to the  value  "value" (typically  0.0, less  often
  1.0). With the value of zpad equal  to the (half) the box size L the
  value of the local variable "border" is 1.
*/
void bgy3d_boundary_set (const State *BHD, Vec g, real value)
{
  const Boundary vol = make_boundary (BHD->PD);

  /* Loop over local portion of grid: */
  {
    /* Get local portion of the grid */
    int x[3], n[3];
    DAGetCorners (BHD->da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

    PetscScalar ***g_;
    DAVecGetArray (BHD->da, g, &g_);

    /*
      The condition in the body of the loop is the opposite of:

      border < i < N[0] - border
      border < j < N[1] - border
      border < k < N[2] - border
    */
    for (int k = x[2]; k < x[2] + n[2]; k++)
      for (int j = x[1]; j < x[1] + n[1]; j++)
        for (int i = x[0]; i < x[0] + n[0]; i++)
          if (unlikely (!inside_boundary (&vol, i, j, k)))
            g_[k][j][i] = value;

    DAVecRestoreArray (BHD->da, g, &g_);
  }
}

