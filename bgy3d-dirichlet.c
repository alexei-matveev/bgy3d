/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O_solutes.c,v 1.3 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

#include <stdbool.h>
#include "bgy3d.h"
#include "bgy3d-vec.h"          /* bgy3d_da_ref() */
#include "bgy3d-dirichlet.h"

// #define MATRIX_FREE

typedef struct Boundary
{
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
 *
 * Convenience   types   and   functions   to   create   matrix   free
 * implementations of linear operators:
 *
 */

/* Functions that take a Mat and a Vec and update another Vec: */
typedef PetscErrorCode (*Operation)(Mat A, Vec x, Vec y);

/* Functions that take a Mat and destroy them: */
typedef PetscErrorCode (*Destructor)(Mat A);


/* Untyped  convinience wrapper,  not all  matrices have  a meaningful
   context: */
static void* context (Mat A)
{
  void *ctx;
  MatShellGetContext (A, &ctx);
  return ctx;
}


/* Creates a matix  shell, but does not associate  any operations with
   it: */
static Mat mat_create_shell (int n, int N, void *ctx)
{
  Mat A;
  MatCreateShell (PETSC_COMM_WORLD, n, n, N, N, ctx, &A);
  return A;
}


static Mat mat_create (int n, int N, void *ctx,
                       Operation mat_mult,
                       Destructor mat_destroy)
{
  /* Create  matrix shell  with  proper dimensions  and associate  the
     context with it: */
  Mat A = mat_create_shell (n, N, ctx);

  /* Set matrix operations: */
  MatShellSetOperation (A, MATOP_MULT,
                        (void (*)(void)) mat_mult);

  MatShellSetOperation (A, MATOP_DESTROY,
                        (void (*)(void)) mat_destroy);

  if (verbosity > 0)
    printf ("mat_create(%p)\n", A);
  return A;
}


/*
 *
 * Two   implementations  of   the  linear   operator   for  Dirichlet
 * problem. The linear operator  that is Laplacian "almost" everywhere
 * except on the boundary.
 *
 */

/*
  The  alternative works  on a  single CPU  too. Again  the regression
  tests  "fail"  with 3-4  digits  the  same  as before.   The  bigger
  difference  appears  to be  that  the matrix  in  AIJ  format has  a
  dedicated  solver for  linear equations.   The  algorithm KSPSolve()
  uses with matrix free implementation  appears to be slower and fails
  for more than one CPU.
*/
#ifndef MATRIX_FREE
/*
  Create and initialize Laplace matrix with appropriate stencil.  This
  may be  implemented using matrix  free facilities of  Petsc, however
  the  one  that  was  tried  is  slower.  The  interface  has  to  be
  consistent with the other impl (see #else):
*/
static Mat lap_mat_create (const DA da, const real h[3],
                           const Boundary *const vol)
{
  const PetscScalar one = 1.0;

  if (verbosity > 0)
    PetscPrintf (PETSC_COMM_WORLD, "Assembling Matrix...");

  /* Create Matrix with appropriate non-zero structure */
  Mat M;
  DAGetMatrix (da, MATMPIAIJ, &M);

  MatZeroEntries (M);

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
              MatSetValuesStencil (M, 1, &row, 1, &row, &one, ADD_VALUES);
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

                  MatSetValuesStencil (M, 1, &row, 3, col, v, ADD_VALUES);
                }
            }
        }

  MatAssemblyBegin (M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd (M, MAT_FINAL_ASSEMBLY);

  if (verbosity > 0)
    PetscPrintf (PETSC_COMM_WORLD, "done.\n");

  return M;
}

#else

/*
  For a  plain Laplacian only the  grid dimensions and  mesh sizes are
  necessary. For a Laplacian with boundary also the third field is set
  and used:
*/
typedef struct Operator
{
  DA da;                        /* array descriptor */
  real h[3];                    /* grid spacing */
  Boundary vol;                 /* only for boundary problems */
} Operator;


static Operator* op_create (const DA da, const real h[3],
                            const Boundary *vol)
{
  /* Make sure it is freed in op_destroy(): */
  Operator *op = malloc (sizeof *op);

  /*
    That is all we need  for the Laplace operator, dimensions and mesh
    sizes.   Copy contents,  the  there are  no  guarantees about  the
    life-time of the object pointed to:
  */
  for (int i = 0; i < 3; i ++)
    op->h[i] = h[i];

  /* We will  DADestroy() it  in op_destroy(), but  it came  from user
     space, so increment the refcount: */
  op->da = bgy3d_da_ref (da);

  /* This is used to copy  over the boundary values: */
  if (vol)
    op->vol = *vol;

  return op;
}


static void op_destroy (Operator *op)
{
  DADestroy (op->da);           /* destroys or decrements refcount */

  free (op);
}


static PetscErrorCode mat_destroy_op (Mat A)
{
  if (verbosity > 0)
    printf ("mat_destroy_op(%p)\n", A);

  /* Only Operators are accepted here: */
  Operator *op = context (A);

  op_destroy (op);

  return 0;
}


/* Does y = A * x = Δx. Laplace operator by finite differences: */
static PetscErrorCode mat_mult_op_lap (Mat A, Vec x, Vec y)
{
  /* Only  matrices  constructed   by  mat_create_lap()  are  accepted
     here: */
  const Operator *lap = context (A);

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
static PetscErrorCode mat_mult_op_bnd (Mat A, Vec x, Vec y)
{
  mat_mult_op_lap (A, x, y);    /* y := Δx */

  Operator *lap = context (A);

  /* Untill  now we compute  plain-old Δx,  now the  boundary specific
     staff. Copy the values at the boundary from x to y as is: */
  copy_boundary (lap->da, &lap->vol, x, y);

  return 0;
}


/*
  To create a matrix one needs  the local and total size of the matrix
  (section). This  is how  to get it  from the array  descriptor.  See
  also msizes().
*/
static void asizes (const DA da, int *n3, int *N3)
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

  *n3 = n[2] * n[1] * n[0];
  *N3 = N[2] * N[1] * N[0];
}


/* Creates Laplacian  matrix. The interface has to  be consistent with
   the other impl (see #ifndef): */
static Mat lap_mat_create (const DA da, const real h[3],
                           const Boundary *vol)
{
  /* Allocates storage for a Operator struct. Make sure to it is freed
     in mat_destroy(): */
  Operator *op = op_create (da, h, vol);

  /* Get the shape of the future matrix: */
  int n, N;
  asizes (da, &n, &N);

  Mat A;
  /* Build either a true Laplacian operator or the one adapted for the
     boundary problem: */
  if (vol)
    A = mat_create (n, N, op, mat_mult_op_bnd, mat_destroy_op);
  else
    A = mat_create (n, N, op, mat_mult_op_lap, mat_destroy_op);

  return A;
}
#endif  /* ifndef MATRIX_FREE */

/*
 *
 * Matrix-free implementation of the matrix inverse.
 *
 */

/* KSP  stays for  Krylov Sub-Space,  used to  solve system  of linear
   equations: */
static KSP ksp_create (Mat M)
{
  PC pc;
  KSP ksp;

  /* Create ksp environment */
  KSPCreate (PETSC_COMM_WORLD, &ksp);
  KSPGetPC (ksp, &pc);

  /* FIXME: literal tolerances here: */
  KSPSetTolerances (ksp, 1.0e-4, 1.0e-4, 1.0e+5, 1000);

  /* Set Matrix */
  //KSPSetOperators (*ksp, M, M, SAME_NONZERO_PATTERN);
  KSPSetOperators (ksp, M, M, SAME_PRECONDITIONER);

  /* Set preconditioner */
  PCSetType (pc, PCBJACOBI);

  /* This is the place which is  supposed to tell KSP solver to re-use
     the supplied vector as initial guess: */
  KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);

  /* runtime options will override default parameters */
  //KSPSetFromOptions(BHD->ksp);

  return ksp;
}


/*
  To create a matrix one needs  the local and total size of the matrix
  (section).  This  is how  to get it  from another matrix.   See also
  asizes()
*/
static void msizes (const Mat A, int *n3, int *N3)
{
  int m, n;
  MatGetLocalSize (A, &m, &n);
  assert (m == n);

  int M, N;
  MatGetSize (A, &M, &N);
  assert (M == N);

  *n3 = n;
  *N3 = N;
}


/* Does y = A x where A  = B^-1 and was constructed as A = mat_inverse
   (B): */
static PetscErrorCode mat_mult_inv (Mat A, Vec x, Vec y)
{
  KSP ksp = context (A);

  KSPSolve (ksp, x, y);

  /* If you think you need to know how many iteration were required to
     solve the equation do this: */
  if (verbosity > 0)
    {
      int iter;
      KSPGetIterationNumber (ksp, &iter);
      PetscPrintf (PETSC_COMM_WORLD, "ksp(%2d) ", iter);
    }

  return 0;
}


static PetscErrorCode mat_destroy_inv (Mat A)
{
  if (verbosity > 0)
    printf ("mat_destroy_inv(%p)\n", A);

  KSP ksp = context (A);

  KSPDestroy (ksp);

  return 0;
}


static Mat mat_inverse (Mat A)
{
  /* An  inverse operator  needs to  solve linear  equations  with the
     original matrix.  Presumable KSP object will hold a reference: */
  KSP ksp = ksp_create (A);

  /* Get the shape of the future matrix: */
  int n, N;
  msizes (A, &n, &N);

  return mat_create (n, N, ksp, mat_mult_inv, mat_destroy_inv);
}


typedef struct Dirichlet
{
  DA da;          /* Array descriptor */
  Mat A;          /* Inverse of the "almost" Laplacian */
  Vec b;          /* Work vector to hold projection on the boundary */
  Boundary vol;   /* Boundary definition */
} Dirichlet;


static PetscErrorCode mat_destroy_dir (Mat A)
{
  if (verbosity > 0)
    printf ("mat_destroy_dir(%p)\n", A);

  Dirichlet *op = context (A);

  DADestroy (op->da);
  MatDestroy (op->A);
  VecDestroy (op->b);

  free (op);

  return 0;
}


static PetscErrorCode mat_mult_dir (Mat L, Vec v, Vec x)
{
  Dirichlet *op = context (L);

  /*
    Get boundary b of v, the rest  of b is set to zero. Together it is
    a linear, albeit not invertible projection operation:

    b := P * v
  */
  VecSet (op->b, 0.0);
  copy_boundary (op->da, &op->vol, v, op->b);

  /*
    Solve Laplace  equation, update  x iteratively.  Ideally  from the
    state of the previous iteration:
            -1
    x := KSP   *  b
  */
  MatMult (op->A, op->b, x);

  /* If you preserve the value of  x until the next call the iterative
     solver will re-use it as initial approximation for the next x. */

  return 0;
}


/* Assemble Laplacian matrix and create KSP environment: */
Mat bgy3d_dirichlet_create (const DA da, const ProblemData *PD)
{
  Dirichlet *op = malloc (sizeof *op);

  const Boundary vol = make_boundary (PD);

  /* I created you ... */
  Mat B = lap_mat_create (da, PD->h, &vol);

  op->da = bgy3d_da_ref (da);    /* DADestroy() it! */
  op->A = mat_inverse (B);       /* MatDestroy() it! */
  op->b = bgy3d_vec_create (da); /* VecDestroy() it! */
  op->vol = vol;

  /* ...   so   I  will   destroy  you  too   (Mat  A   holds  another
     reference): */
  MatDestroy (B);

  /* Get the shape of the future matrix: */
  int n, N;
  msizes (op->A, &n, &N);

  return mat_create (n, N, op, mat_mult_dir, mat_destroy_dir);
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
void bgy3d_impose_laplace_boundary (const State *BHD, Vec v, Vec x)
{
  /*
    Get boundary b of v, the rest  of b is set to zero. Together it is
    a linear, albeit not invertible projection operation:

    b := P * v

    and solve  Laplace equation,  update x iteratively.   Ideally from
    the state of the previous iteration:

            -1
    x := KSP   * b

    The Mat dirichlet_mat implements these two linear operations:
  */
  MatMult (BHD->dirichlet_mat, v, x);

  /*
    Now subtract  solution from  v.  The whole  is in effect  a linear
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
  Boundary vol = make_boundary (BHD->PD);
  set_boundary (BHD->da, &vol, v, 0.0);

  /* If you preserve the value of  x until the next call the iterative
     solver will re-use it as initial approximation for the next x. */
}

/* Laplace matrix, no boundary: */
Mat bgy3d_lap_mat_create (const DA da, const real h[3])
{
  return lap_mat_create (da, h, NULL);
}

