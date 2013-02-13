/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/*
 * Provides  a   Petsc  Mat  interface   to  the  FFTW   package.  See
 * bgy3d_fft_test() for an example use.
 */

#include "bgy3d.h"              /* KFREQ(), M_PI */
#include <assert.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include "petscda.h"            /* DA, Vec */
#include "bgy3d-vec.h"          /* bgy3d_da_ref() */
#include "bgy3d-fftw.h"         /* Common interface for two impls */
#include <complex.h>            /* after fftw.h */

typedef struct {
  /* Array  descriptors for real  and complex  vectors that  share the
     distribution pattern with FFTW-MPI: */
  DA da, dc;

  /* Two plans for doubl -> cmplx and cmplx -> doubl: */
  fftw_plan fw, bw;

  /*
    Storage for a PADDED array  of reals and complex numbers. The last
    dimension of the double array is  2 (N/2 + 1), that of the complex
    array just N/2 + 1:
  */
  double *doubl;
  fftw_complex *cmplx;
} FFT;

/*
  Always  use  real  array   descriptor  to  get  global  grid  shape.
  Dimensions  of complex  arrays  maybe smaller.  Shape  of the  local
  section is yet another thing.
*/
static void shape (const FFT *fft, int *NI, int *NJ, int *NK)
{
  int dim;
  DAGetInfo (fft->da, &dim, NI, NJ, NK,
             NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  assert (dim == 3);
}

static const int debug = 0;

/* doubl := Vec, for forward FFT */
static void unpack_real (FFT *fft, Vec g, double *restrict doubl)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (fft->da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the grid: */
  shape (fft, &NI, &NJ, &NK);

  /* Padded (physically the last) dimension. */
  const int nip = 2 * (NI / 2 + 1);
  if (debug)
    printf ("unpack_real: shape = %d %d %d\n", nk, nj, nip);
  assert (ni < nip);

  /* The  view of  local FFT  padded storage  as a  3d array  with the
     proper highest dimension: */
  double (*const view)[nk][nj][nip] =(double (*)[nk][nj][nip]) doubl;

  /* loop over local portion of grid */
  {
    double ***g_;
    DAVecGetArray (fft->da, g, &g_);

    /* NOTE: array pointed at by view is padded: */
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          (*view)[k][j][i] = g_[k0 + k][j0 + j][i0 + i];

    DAVecRestoreArray (fft->da, g, &g_);
  }
}


/* Vec := doubl, for inverse FFT */
static void pack_real (FFT *fft, Vec g, const double *restrict doubl)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (fft->da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the grid: */
  shape (fft, &NI, &NJ, &NK);

  /* Padded (physically the last) dimension. */
  const int nip = 2 * (NI / 2 + 1);
  if (debug)
    printf ("pack_real: shape = %d %d %d\n", nk, nj, nip);
  assert (ni < nip);

  /* The  view of  local FFT  padded storage  as a  3d array  with the
     proper highest dimension: */
  double (*const view)[nk][nj][nip] =(double (*)[nk][nj][nip]) doubl;

  /* loop over local portion of grid */
  {
    double ***g_;
    DAVecGetArray (fft->da, g, &g_);

    /* NOTE: array pointed at by view is padded: */
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          g_[k0 + k][j0 + j][i0 + i] = (*view)[k][j][i];

    DAVecRestoreArray (fft->da, g, &g_);
  }
}


/* Vec := cmplx, for forward FFT */
static void pack_cmplx (FFT *fft, Vec g, /* const */ fftw_complex *cmplx)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (fft->dc, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the grid: */
  shape (fft, &NI, &NJ, &NK);

  /* Padded (physically the last) dimension: */
  const int nip = NI / 2 + 1;
  if (debug)
    printf ("pack_cmplx: shape = %d %d %d\n", nk, nj, nip);

  /*
    Otherwise  we  cannot  convert  between complex  and  half-complex
    locally (this is probably also the  reason why there is no real to
    half-complex transforms in FFTW3 MPI):
  */
  assert (ni == nip);
  assert (i0 == 0);

  /* The view  of local  FFT complex  storage as a  3d array  with the
     proper highest dimension: */
  complex (*const view)[nk][nj][nip] = (complex (*)[nk][nj][nip]) cmplx;

  /* loop over local portion of grid */
  {
    complex ***g_;

    DAVecGetArray (fft->dc, g, &g_);

    /* loop over local portion of grid */
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < nip; i++)
          g_[k0 + k][j0 + j][i0 + i] = (*view)[k][j][i];

    DAVecRestoreArray (fft->dc, g, &g_);
  }
}


/* cmplx := Vec, for inverse FFT */
static void unpack_cmplx (FFT *fft, Vec g, fftw_complex *cmplx)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (fft->dc, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the grid: */
  shape (fft, &NI, &NJ, &NK);

  /* Padded (physically the last) dimension: */
  const int nip = NI / 2 + 1;
  if (debug)
    printf ("unpack_cmplx: shape = %d %d %d\n", nk, nj, nip);
  assert (ni < 2 * nip);

  /*
    Otherwise  we  cannot  convert  between complex  and  half-complex
    locally (this is probably also the  reason why there is no real to
    half-complex transforms in FFTW3 MPI):
  */
  assert (ni == nip);
  assert (i0 == 0);

  /* The view  of local  FFT complex  storage as a  3d array  with the
     proper highest dimension: */
  complex (*const view)[nk][nj][nip] = (complex (*)[nk][nj][nip]) cmplx;

  /* loop over local portion of grid */
  {
    complex ***g_;

    DAVecGetArray (fft->dc, g, &g_);

    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < nip; i++)
          (*view)[k][j][i] = g_[k0 + k][j0 + j][i0 + i];

    DAVecRestoreArray (fft->dc, g, &g_);
  }
}


static FFT* context (Mat A)
{
  FFT *fft;
  MatShellGetContext (A, (void**) &fft);
  return fft;
}


/* Does y = A * x. Forward FFT with the interface for use by Petsc: */
static PetscErrorCode mat_mult_fft (Mat A, Vec x, Vec y)
{
  /* Only  matrices  constructed   by  mat_create_fft()  are  accepted
     here: */
  FFT *fft = context (A);

  /* Fill real array with real data from x: */
  unpack_real (fft, x, fft->doubl);

  /* forward fft */
  fftw_execute (fft->fw);

  /* Pack complex output into complex Vec y: */
  pack_cmplx (fft, y, fft->cmplx);

  return 0;
}


/* Does y = A^T * x. The inverse FFT. */
static PetscErrorCode mat_mult_transpose_fft (Mat A, Vec x, Vec y)
{
  /* Only matrices constructed  by bgy3d_fft_mat_create() are accepted
     here: */
  FFT *fft = context (A);

  /* Fill complex array with halfcomplex data from x: */
  unpack_cmplx (fft, x, fft->cmplx);

  /* inverse fft */
  fftw_execute (fft->bw);

  /* Pack real output into Vec y: */
  pack_real (fft, y, fft->doubl);

  return 0;
}


static PetscErrorCode mat_destroy_fft (Mat A)
{
  /* Only matrices constructed  by bgy3d_fft_mat_create() are accepted
     here: */
  FFT *fft = context (A);

  /*
    Since  bgy3d_fft_mat_create()  returns  DA  for real  and  complex
    vectors  to the caller,  it is  his/her responsibility  to destroy
    them after  use.  We  do this  here too. To  not let  him/her with
    dangling pointers we should have incremented the reference counter
    on these two objects in bgy3d_fft_mat_create().
  */
  DADestroy (fft->da);
  DADestroy (fft->dc);

  fftw_destroy_plan (fft->fw);
  fftw_destroy_plan (fft->bw);

  fftw_free (fft->doubl);
  fftw_free (fft->cmplx);
  free (fft);

  return 0;
}


/*
  Create internals of  the matrix.  There is a  mind bending caveat to
  keep in  mind. When the dimensions  N[3] = {5, 7,  11} then whenever
  you  think  of memory  layout  think  of  a column  major  (fortran)
  f_layout(5, 7, 11) or  c_layout[11][7][5].  Petsc advertizes the use
  of  the call sequence  suggesting that  the array  NI x  NJ x  NK is
  stored in the column major order with the stride-1 dimension NI:

    DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);
    DAVecGetArray (da, vec, &layout);

    int ijk = 0;
    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          buf[ijk++] = layout[k][j][i];

  So  to possibly decrease  the impedance  mismatch between  Petsc and
  FFTW we  declare the first  array element N[0]  to the the  one with
  unit stride.  Then, of course N[2]  is the major  dimension with the
  largest stride. That is why the order of N[2], N[1], and N[0] in the
  calls to FFTW routines is reversed.

  A consequence of that is that the work sharing is along N[2] as this
  is the "leading"  or "major" dimension in the  physical layout. This
  cannot be changed, given the above  choice, as FFTW MPI has only one
  way to share the work. For  the above example 11 major columns would
  have  been distributed  among the  workers,  say as  6 +  5 for  two
  workers.

  Petsc, on the  other hand is flexible with  data distribution layout
  and will  adapt its  array descriptors to  distribute the  data over
  N[2]  too.  N[2] given  the "column  major" interpretation  of Petsc
  arrays is also the "major" dimension so that both FFTW and Petsc are
  consistent in that respect.

  And,  finally, there  is  a padding  issue.   The stride-1  physical
  dimension of the FFTW arrays deserves a special care. To ever pass a
  real array real[11][7][5] that  corresponds to our interpretation of
  N[3]  =  {5,   7,  11}  you  will  need  a   storage  layed  out  as
  buf[11][7][6].   Note the  padding in  the rightmost  position.  The
  padded dimension is computed as 6 = 2 (5/2 + 1). Or, in general NP =
  2  * (N/2 +  1).  The  corresponding complex  arrays will  have NP/2
  elements which is the main reason for the padding, actually.
*/
static int fftw_mpi_init_called = 0;
void bgy3d_fft_mat_create (const int N[3], Mat *A, DA *da, DA *dc)
{
  /* FIXME: find a better place: */
  if (!fftw_mpi_init_called)
    {
      fftw_mpi_init ();          /* required */
      atexit (fftw_mpi_cleanup); /* required */
      fftw_mpi_init_called = 1;
    }

  /* Allocates storage for an FFT struct: */
  FFT *fft = malloc (sizeof *fft);

  /* Get number of processes */
  int np, id;
  MPI_Comm_size (PETSC_COMM_WORLD, &np);
  MPI_Comm_rank (PETSC_COMM_WORLD, &id);

  /* FIXME: see bgy3d_fft_init_da () and avoid code duplication: */
  {
    ptrdiff_t alloc_local, local_range, local_start;

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_3d (N[2], N[1], N[0] / 2 + 1,
                                          PETSC_COMM_WORLD,
                                          &local_range,
                                          &local_start);

    /* Scratch arrays for FFT: */
    fft->doubl = fftw_alloc_real (2 * alloc_local);
    fft->cmplx = fftw_alloc_complex (alloc_local);

    /* create plan for out-of-place forward DFT */
    fft->fw = fftw_mpi_plan_dft_r2c_3d (N[2], N[1], N[0],
                                        fft->doubl, fft->cmplx,
                                        PETSC_COMM_WORLD,
                                        FFTW_ESTIMATE);
    assert (fft->fw != NULL);

    /* create plan for out-of-place inverse DFT */
    fft->bw = fftw_mpi_plan_dft_c2r_3d (N[2], N[1], N[0],
                                        fft->cmplx, fft->doubl,
                                        PETSC_COMM_WORLD,
                                        FFTW_ESTIMATE);
    assert (fft->bw != NULL);

    /*
      Create  Petsc  Distributed  Array  according  to  FFTW-MPI  data
      distribution. The  FFTW MPI distributes the  work/data among the
      workers along the leading dimension:
    */
    int l0[1], l1[1], l2[np]; /* sum (l2[:]) == N[0] */

    /* Collect  the  info  about  the  local  shares  of  the  leading
       dimension from all workers: */
    {
      int local_range_ = local_range; /* cast to int */
      int err = MPI_Allgather (&local_range_, 1, MPI_INT, l2, 1, MPI_INT, PETSC_COMM_WORLD);
      assert (err == MPI_SUCCESS);
    }
    l0[0] = N[0];
    l1[0] = N[1];

    if (debug)
      {
        printf ("%d: local_range = %d, local_start = %d\n",
                id, (int) local_range, (int) local_start);

        printf ("%d: alloc_local = %d\n", id, (int) alloc_local);

        for (int p = 0; p < np; p++)
          printf ("%d: l0[%d] = %d\n", id, p, l2[p]);
      }

    const PetscInt stencil_width = 1;

    /*
      Petsc  appears to  use different  convention similar  to fortran
      column major. So that the leading dimension along which the work
      is shared appears  last here. As of Petsc  3 the library refuses
      (loudly) to handle  zero range for any worker.   This limits the
      maximum number of workers by N[2]:
    */
    DACreate3d (PETSC_COMM_WORLD,
                DA_XYZPERIODIC,    /* was DA_NONPERIODIC */
                DA_STENCIL_STAR,
                N[0], N[1], N[2],
                1, 1, np,
                1, stencil_width,
                l0, l1, l2,
                &fft->da);

    /* For complex vectors: */
    int l0half[1] = {N[0] / 2 + 1};
    DACreate3d (PETSC_COMM_WORLD,
                DA_XYZPERIODIC,    /* was DA_NONPERIODIC */
                DA_STENCIL_STAR,
                N[0] / 2 + 1, N[1], N[2], /* ~half as many, but ... */
                1, 1, np,
                2, stencil_width, /* ... with two DOF */
                l0half, l1, l2,
                &fft->dc);
  }

  /* Get local dimensions: */
  int x[3], n[3];
  DAGetCorners (fft->da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  if (debug)
    {
      for (int i = 0; i < 3; i++)
        printf ("%d: n[%d] = %d, x[%d] = %d\n",
                id, i, n[i], i, x[i]);
    }

  /* This is not true anymore as FFTW MPI needs some padding: */
  /* assert (n[0] * n[1] * n[2] == alloc_local); */

  /*
    Create the matrix itself:

    First, set  total and local size  of the matrix  (section). I gues
    this  only works  because the  work/storage is  divided  along one
    dimension:
  */
  {
    const int M0 = N[2] * N[1] * (N[0] / 2 + 1) * 2;
    const int M1 = N[2] * N[1] * N[0];

    const int m0 = n[2] * n[1] * (n[0] / 2 + 1) * 2;
    const int m1 = n[2] * n[1] * n[0];

    if (debug)
      {
        printf ("%d: matrix dim %d x %d (local), %d x %d (global)\n",
                id, m0, m1, M0, M1);
      }
    /* Also puts internals (FFT struct) into the matrix: */
    MatCreateShell (PETSC_COMM_WORLD, m0, m1, M0, M1, (void*) fft, A);
  }

  /* Set matrix operations: */
  MatShellSetOperation (*A, MATOP_MULT,
                        (void (*)(void)) mat_mult_fft);

  MatShellSetOperation (*A, MATOP_MULT_TRANSPOSE,
                        (void (*)(void)) mat_mult_transpose_fft);

  MatShellSetOperation (*A, MATOP_DESTROY,
                        (void (*)(void)) mat_destroy_fft);

  /*
    Also return  DA descriptors for  real and complex vectors  so that
    the user can create them.   But first increase the reference count
    on  these  two objects  because  the  FFT  struct refers  to  them
    too. The  user code  is still responsible  to call  DADestroy() on
    both objects.  Otherwise the refcount will never reach zero:
  */
  *da = bgy3d_da_ref (fft->da);
  *dc = bgy3d_da_ref (fft->dc);
}


/*
  Given the FFT transform Vec  Y of some equally-spaced grid sample of
  a real-valued function y(x), compute the values of y(x) at arbitrary
  given  points x[][]  by trigonometric  interpolation. Note  that the
  scale of  x here is such  that the components range  between 0.0 and
  the respective N[] within  the unit cell.  Extrapolation outside the
  cell is periodic by the nature of trigonometric interpolation.

  Note that "upscaling"  to a finer but regular grid  can be done more
  efficiently.
*/
void bgy3d_fft_interp (const Mat A,
                       const Vec Y, /* complex, intent(in) */
                       int np, double x[np][3], /* intent(in) */
                       double y[np])            /* intent(out) */
{
  const FFT *fft = context (A);

  /* Get true grid shape N[3]: */
  int N[3];
  shape (fft, &N[0], &N[1], &N[2]);

  /* Init    real   accumulator.    Imaginary    part   vanishes    by
     construction: */
  for (int p = 0; p < np; p++)
    y[p] = 0.0;

  /* loop over local portion of grid */
  int c[3], n[3];             /* corner and the size of the section */
  DAGetCorners (fft->dc, &c[0], &c[1], &c[2], &n[0], &n[1], &n[2]);

  complex ***Y_;
  DAVecGetArray (fft->dc, Y, &Y_);

  int k[3];
  for (k[2] = c[2]; k[2] < c[2] + n[2]; k[2]++)
    for (k[1] = c[1]; k[1] < c[1] + n[1]; k[1]++)
      for (k[0] = c[0]; k[0] < c[0] + n[0]; k[0]++)
        {
          /* Take negative frequencies where k > N/2: */
          int K[3];
          FOR_DIM
            K[dim] = KFREQ (k[dim], N[dim]);

          /* Note the sequence, 2, 1, 0: */
          const complex yk = Y_[k[2]][k[1]][k[0]];

          for (int p = 0; p < np; p++)
            {
              /* Coordinates  x[p][]  are  fractional, represented  by
                 real numbers: */
              complex wkx = 1.0;
              FOR_DIM
                if ((N[dim] - 2 * K[dim]) % N[dim] != 0) /* k != N - k mod N */
                  wkx *= cexp (2 * M_PI * K[dim] * x[p][dim] / N[dim] * I);
                else
                  wkx *= cos (2 * M_PI * K[dim] * x[p][dim] / N[dim]);

              /* Avoid double counting: */
              if ((N[0] - 2 * K[0]) % N[0] == 0) /* k == N - k mod N */
                wkx /= 2;

              /*
                As  only (roughly)  a half  of complex  amplitudes are
                stored, add two terms for  +K and -K. FIXME: This term
                is  real  by  construction  maybe use  that  to  avoid
                complex arithmetics?
              */
              y[p] += yk * wkx + conj (yk * wkx);
            }
        }
  DAVecRestoreArray (fft->dc, Y, &Y_);

  /* Each  worker summed  only over  its own  range of  K[],  sum over
     workers: */
  bgy3d_comm_allreduce (np, y);

  /*
    FIXME:  *do*  normalize  as  for  interpolation.   This  does  not
    correspond  to the  inverse FFT  as implemented  by FFTW  which is
    unnormalized!
  */
  for (int p = 0; p < np; p++)
    y[p] /= (N[0] * N[1] * N[2]);
}
