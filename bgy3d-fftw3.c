/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */

/*
 * Provides  a   Petsc  Mat  interface   to  the  FFTW   package.  See
 * bgy3d_test_fft() for an example use.
 */

// #include <unistd.h>             /* getpid(), remove later */
#include <assert.h>
// #include <complex.h>          /* makes setting re and im cumbersome */
#include <fftw3.h>
#include <fftw3-mpi.h>
#include "petscda.h"            /* DA, Vec */
#include "bgy3d-fftw3.h"

typedef struct {
  /* Array  descriptor  that  shares  the  distribution  pattern  with
     FFTW-MPI: */
  DA da;

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


/* doubl := Vec: */
static void unpack_real (DA da, Vec g, double *restrict doubl)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the array: */
  DAGetInfo (da, NULL, &NI, &NJ, &NK,
             NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  /* Padded (physically the last) dimension. */
  const int nip = 2 * (NI / 2 + 1);
  /* printf ("unpack_real: shape = %d %d %d\n", nk, nj, nip); */
  assert (ni < nip);

  /* The  view of  local FFT  padded storage  as a  3d array  with the
     proper highest dimension: */
  double (*const view)[nk][nj][nip] =(double (*)[nk][nj][nip]) doubl;

  /* loop over local portion of grid */
  {
    PetscScalar ***g_vec;
    DAVecGetArray (da, g, &g_vec);

    /* NOTE: array pointed at by view is padded: */
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          (*view)[k][j][i] = g_vec[k0 + k][j0 + j][i0 + i];

    DAVecRestoreArray (da, g, &g_vec);
  }
}


/* Vec := doubl: */
static void pack_real (DA da, Vec g, const double *restrict doubl)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the array: */
  DAGetInfo (da, NULL, &NI, &NJ, &NK,
             NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  /* Padded (physically the last) dimension. */
  const int nip = 2 * (NI / 2 + 1);
  /* printf ("pack_real: shape = %d %d %d\n", nk, nj, nip); */
  assert (ni < nip);

  /* The  view of  local FFT  padded storage  as a  3d array  with the
     proper highest dimension: */
  double (*const view)[nk][nj][nip] =(double (*)[nk][nj][nip]) doubl;

  /* loop over local portion of grid */
  {
    PetscScalar ***g_vec;
    DAVecGetArray (da, g, &g_vec);

    /* NOTE: array pointed at by view is padded: */
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        for (int i = 0; i < ni; i++)
          g_vec[k0 + k][j0 + j][i0 + i] = (*view)[k][j][i];

    DAVecRestoreArray (da, g, &g_vec);
  }
}


/* Vec := cmplx: */
static void pack_cmplx (DA da, Vec g, /* const */ fftw_complex *cmplx)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the array: */
  DAGetInfo (da, NULL, &NI, &NJ, &NK,
             NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  /*
    Otherwise  we  cannot  convert  between complex  and  half-complex
    locally (this is probably also the  reason why there is no real to
    half-complex transforms in FFTW3 MPI):
  */
  assert (ni == NI);
  assert (i0 == 0);

  /* Padded (physically the last) dimension: */
  const int nip = NI / 2 + 1;
  /* printf ("pack_cmplx: shape = %d %d %d\n", nk, nj, nip); */

  /* The view  of local  FFT complex  storage as a  3d array  with the
     proper highest dimension: */
  fftw_complex (*const view)[nk][nj][nip] = (fftw_complex (*)[nk][nj][nip]) cmplx;

  /* loop over local portion of grid */
  {
    PetscScalar ***g_vec;
    DAVecGetArray (da, g, &g_vec);

    /* loop over local portion of grid */
    /* Attention: order of indices is not variable */
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        {
          /*
            Halfcomplex format: the  real array of size n  is built of
            (non-vanishing)  real and imaginary  parts of  the complex
            array of size  n / 2 +  1. Here is the layout  of the real
            array:

            r, r , ..., r   , i         , ..., i , i
             0  1        n/2   (n+1)/2-1        2   1

            At  positions 0  <= i  <= NI/2 in
            doubl[] one puts the real parts of cmplx[i]:
          */
          for (int i = 0; i <= NI / 2; i++)
            g_vec[k0 + k][j0 + j][i0 + i] = (*view)[k][j][i][0];

          /*
            At  positions  NI/2 <  i  < NI  in  doubl[]  one puts  the
            imaginary parts of cmplx[??]:
          */
          for (int i = NI - 1; i > NI / 2; i--)
            g_vec[k0 + k][j0 + j][i0 + i] = (*view)[k][j][NI - i][1];
        }
    DAVecRestoreArray (da, g, &g_vec);
  }
}


/* cmplx := Vec: */
static void unpack_cmplx (DA da, Vec g, fftw_complex *cmplx)
{
  int i0, j0, k0, ni, nj, nk, NI, NJ, NK;

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* Get the official (full) dimension of the array: */
  DAGetInfo (da, NULL, &NI, &NJ, &NK,
             NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  /*
    Otherwise  we  cannot  convert  between complex  and  half-complex
    locally (this is probably also the  reason why there is no real to
    half-complex transforms in FFTW3 MPI):
  */
  assert (ni == NI);
  assert (i0 == 0);

  /* Padded (physically the last) dimension: */
  const int nip = NI / 2 + 1;
  /* printf ("unpack_cmplx: shape = %d %d %d\n", nk, nj, nip); */

  /* The view  of local  FFT complex  storage as a  3d array  with the
     proper highest dimension: */
  fftw_complex (*const view)[nk][nj][nip] = (fftw_complex (*)[nk][nj][nip]) cmplx;

  /* loop over local portion of grid */
  {
    PetscScalar ***g_vec;
    DAVecGetArray (da, g, &g_vec);

    for (int k = 0; k < nk; k++)
      for (int j = 0; j < nj; j++)
        {
          /*
            Halfcomplex format: see pack_cmplx().
          */
          for (int i = 0; i <= NI / 2; i++)
            (*view)[k][j][i][0] = g_vec[k0 + k][j0 + j][i0 + i];

          for (int i = NI - 1; i > NI / 2; i--)
            (*view)[k][j][NI - i][1] = g_vec[k0 + k][j0 + j][i0 + i];

          /*
            Frequency  zero (i ==  0) and  n /  2 for  even n  has no
            imaginary part and are not stored in halfcomplex format:
          */
          (*view)[k][j][0][1] = 0.0;

          if (NI % 2 == 0)
            (*view)[k][j][NI / 2][1] = 0.0;
        }
    DAVecRestoreArray (da, g, &g_vec);
  }
}


/* Does y = A * x. Forward FFT with the interface for use by Petsc: */
static PetscErrorCode mat_mult_fft (Mat A, Vec x, Vec y)
{
  /* Only  matrices  constructed   by  mat_create_fft()  are  accepted
     here: */
  FFT *fft;
  MatShellGetContext (A, (void**) &fft);

  /* Fill real array with real data from x: */
  unpack_real (fft->da, x, fft->doubl);

  /* forward fft */
  fftw_execute (fft->fw);

  /* Pack complex output into real y in halfcomplex format: */
  pack_cmplx (fft->da, y, fft->cmplx);

  return 0;
}


/* Does y = A^T * x. The inverse FFT. */
PetscErrorCode mat_mult_transpose_fft (Mat A, Vec x, Vec y)
{
  /* Only  matrices  constructed   by  mat_create_fft()  are  accepted
     here: */
  FFT *fft;
  MatShellGetContext (A, (void**) &fft);

  /* Fill complex array with halfcomplex data from x: */
  unpack_cmplx (fft->da, x, fft->cmplx);

  /* inverse fft */
  fftw_execute (fft->bw);

  /* Pack real output into Vec y: */
  pack_real (fft->da, y, fft->doubl);

  return 0;
}


PetscErrorCode mat_destroy_fft (Mat A)
{
  /* Only  matrices  constructed   by  mat_create_fft()  are  accepted
     here: */
  FFT *fft;
  MatShellGetContext (A, (void**) &fft);

  /* Abstraction  leaking  here:  mat_create_fft()  return DA  to  the
     caller, what if he continues to use it? */
  DADestroy(fft->da);

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
PetscErrorCode mat_create_fft (const int N[3], Mat *A, DA *da)
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
    ptrdiff_t l0[1], l1[1], l2[np]; /* sum (l2[:]) == N[0] */

    /* Collect  the  info  about  the  local  shares  of  the  leading
       dimension from all workers: */
    {
      assert (sizeof(ptrdiff_t) == sizeof(long));
      int err = MPI_Allgather (&local_range, 1, MPI_LONG, l2, 1, MPI_LONG, PETSC_COMM_WORLD);
      assert (err == MPI_SUCCESS);
    }
    l0[0] = N[0];
    l1[0] = N[1];

    if (0)
      {
        printf ("%d: local_range = %d, local_start = %d\n",
                id, local_range, local_start);

        printf ("%d: alloc_local = %d\n", id, alloc_local);

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
  }

  /* Get local dimensions: */
  int x[3], n[3];
  DAGetCorners (fft->da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  if (0)
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
    this  only works  because the  dork/storage is  divided  along one
    dimension:
  */
  {
    int M = 1, m = 1;
    for (int i = 0; i < 3; i++)
      {
        M *= N[i];
        m *= n[i];
      }

    /* Also puts internals into the matrix. FIXME: local sizes? */
    MatCreateShell (PETSC_COMM_WORLD, m, m, M, M, (void*) fft, A);
  }

  /* Set matrix operations: */
  MatShellSetOperation (*A, MATOP_MULT,
                        (void (*)(void)) mat_mult_fft);

  MatShellSetOperation (*A, MATOP_MULT_TRANSPOSE,
                        (void (*)(void)) mat_mult_transpose_fft);

  MatShellSetOperation (*A, MATOP_DESTROY,
                        (void (*)(void)) mat_destroy_fft);

  /* Also return  DA descriptor so  that the user can  create suitable
     vectors: */
  *da = fft->da;

  return 0;
}

double bgy3d_test_fft (int m, int n, int k)
{
  const int N[3] = {m, n, k};
  const int NNN = N[0] * N[1] * N[2];
  Mat A;
  DA da;

  /* pid_t pid = getpid(); */
  /* printf ("PID=%d\n", pid); */
  /* sleep (20); */

  mat_create_fft (N, &A, &da);

  Vec x, y, z;

  DACreateGlobalVector (da, &x);
  DACreateGlobalVector (da, &y);
  DACreateGlobalVector (da, &z);

  VecSet (x, 1.0);
  VecSet (y, 99.0);
  VecSet (z, 999.0);

  /* This corresponds to direct FFT, y = fft(x): */
  MatMult (A, x, y);

  /* VecView (y, PETSC_VIEWER_STDOUT_WORLD); */

  /* This corresponds to inverse FFT, z = fft(y): */
  MatMultTranspose (A, y, z);

  /*
    FIXME: The matrix should probably  be made orthogonal sot that the
    intuitive relation A^T A  * x == x holds. At the  moment A^T * A =
    grid size:
  */
  VecAXPY (z, -NNN, x);

  /* VecView (z, PETSC_VIEWER_STDOUT_WORLD); */
  double norm;
  VecNorm (z, NORM_INFINITY, &norm);

  /* PetscPrintf (PETSC_COMM_WORLD, "norm=%e\n", norm); */

  VecDestroy (x);
  VecDestroy (y);
  VecDestroy (z);

  /* FIXME: Also destroys array descriptor da: */
  MatDestroy (A);

  return norm;
}
