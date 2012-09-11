/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.15 2007-04-23 16:55:13 jager Exp $ */
/*==========================================================*/

#include <assert.h>
#include <fftw_mpi.h>
#include "fft_3d.h"             /* FFT_DATA */
#include "petscda.h"            /* DA, Vec */
#include "bgy3d-fft.h"

extern int verbosity;   /* FIXME: dont want to #include bgy3d.h yet */

/* Initializes  distributes  array descriptors  and  FFTW  plans in  a
   consistent fashion. The data is distributed over in xy-planes (over
   subranges of z-indices).  In the scope of this code  z is the major
   axis, think of array[z][y][x]. */
void bgy3d_fft_init_da (const int N[3],
                        fftwnd_mpi_plan *fw, fftwnd_mpi_plan *bw,
                        DA *da, DA *da_mg)
{
  int x[3], n[3];
  int np, id;
  int local_nx, local_x_start, local_ny, local_y_start, total_local_size;
  PetscInt lx[1], ly[1], *lz;

  /* Initialize parallel stuff: fftw + petsc */
  *fw = fftw3d_mpi_create_plan(PETSC_COMM_WORLD,
                               N[2], N[1], N[0],
                               FFTW_FORWARD, FFTW_ESTIMATE);
  assert (*fw != NULL);

  *bw = fftw3d_mpi_create_plan(PETSC_COMM_WORLD,
                               N[2], N[1], N[0],
                               FFTW_BACKWARD, FFTW_ESTIMATE);

  assert (*bw != NULL);

  fftwnd_mpi_local_sizes(*fw, &local_nx, &local_x_start,
                         &local_ny, &local_y_start, &total_local_size);

  /* Get number of processes */
  MPI_Comm_size (PETSC_COMM_WORLD, &np);
  MPI_Comm_rank (PETSC_COMM_WORLD, &id);

  /* Create Petsc Distributed Array according to fftw data distribution*/
  lz = (PetscInt*) malloc(np*sizeof(*lz));

  MPI_Allgather( &local_nx, 1, MPI_INT, lz, 1, MPI_INT, PETSC_COMM_WORLD);
  ly[0] = N[1];
  lx[0] = N[2];

#if defined(L_BOUNDARY) || defined(L_BOUNDARY_MG)
    const PetscInt stencil_width = 1;
#else
    const PetscInt stencil_width = 0;
#endif

  DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR,
             N[0], N[1], N[2],
             1, 1, np,
             1, stencil_width,
             lx, ly, lz,
             da);

  /* In multigird case, also construct distributed array discriptor of
     half the size. Used only ifdef L_BOUNDARY_MG: */
  if (da_mg) {
      for(int p = 0; p < np; p++)
          lz[p] /= 2;
      lx[0] /= 2;
      ly[0] /= 2;
      DACreate3d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR,
                 N[0] / 2, N[1] / 2, N[2] / 2,
                 1, 1, np,
                 1, stencil_width,
                 lx, ly, lz,
                 da_mg);
  }

  DAGetCorners(*da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  if (verbosity > 2)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Subgrids on processes:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "id %d of %d: %d %d %d\t%d %d %d\tfft: %d %d\n",
                              id, np, x[0], x[1], x[2], n[0], n[1], n[2],
                              local_nx, local_x_start);
      PetscSynchronizedFlush(PETSC_COMM_WORLD);
    }
  assert (n[0] * n[1] * n[2] == total_local_size);

  free(lz);
}

fftw_complex *bgy3d_fft_malloc (DA da)
{
    int x[3], n[3];

    /* Get local portion of the grid */
    DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

    return (fftw_complex*) malloc(n[0] * n[1] * n[2] * sizeof(fftw_complex));
}

void bgy3d_fft_free (fftw_complex *ptr)
{
    free(ptr);
}

static void unpack (DA da, Vec g, fftw_complex *restrict g_fft)
{
    int index, i0, j0, k0, ni, nj, nk;
    PetscScalar ***g_vec;

    /* Get local portion of the grid */
    DAGetCorners(da, &i0, &j0, &k0, &ni, &nj, &nk);

    DAVecGetArray(da, g, &g_vec);

    /* loop over local portion of grid */
    /* Attention: order of indices is not variable */
    index = 0;
    for (int k = k0; k < k0 + nk; k++)
        for (int j = j0; j < j0 + nj; j++)
            for (int i = i0; i < i0 + ni; i++)
                {
                    g_fft[index].re = g_vec[k][j][i];
                    g_fft[index].im = 0;  /* Vec g is real */
                    index++;
                }
    DAVecRestoreArray(da, g, &g_vec);
}

static void pack (DA da, Vec g, const fftw_complex *restrict g_fft)
{
    int index, i0, j0, k0, ni, nj, nk;
    PetscScalar ***g_vec;

    /* Get local portion of the grid */
    DAGetCorners(da, &i0, &j0, &k0, &ni, &nj, &nk);

    DAVecGetArray(da, g, &g_vec);

    /* loop over local portion of grid */
    /* Attention: order of indices is not variable */
    index = 0;
    for (int k = k0; k < k0 + nk; k++)
        for (int j = j0; j < j0 + nj; j++)
            for (int i = i0; i < i0 + ni; i++)
                {
                    /* FIXME: this  is subtle,  we are packing  complex vector
                       into a real array. Imaginary part gets ignored: */
                    g_vec[k][j][i] = g_fft[index].re;
                    index++;
                }
    DAVecRestoreArray(da, g, &g_vec);
}

/*
 * The function  has a feature,  if the output  array g_fft is  NULL a
 * fresh array is allocated with the size derived from the distributed
 * array  description.  This  size  is suffucient  to  hold the  local
 * portion of the array.
 */
fftw_complex *ComputeFFTfromVec_fftw (DA da, fftwnd_mpi_plan fft_plan, Vec g,
				fftw_complex *g_fft, fftw_complex *scratch)
{
    if (g_fft == NULL)
        g_fft = bgy3d_fft_malloc (da);

    /* Real Vec into complex array: */
    unpack (da, g, g_fft);

    /* forward fft */
    fftwnd_mpi( fft_plan, 1, g_fft, scratch, FFTW_NORMAL_ORDER);

    return g_fft;
}


void ComputeVecfromFFT_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
			    fftw_complex *g_fft, fftw_complex *scratch)
{
    assert (g_fft != NULL);

    /* backward fft */
    fftwnd_mpi( fft_plan, 1, g_fft, scratch, FFTW_NORMAL_ORDER);

    /* Pack a  complex vector with (hopefully)  vanishing imaginary part
       into a real Vec: */
    pack (da, g, g_fft);
}

/* y := alpha  * x + beta *  y. FIXME: is there anything  like that in
   FFTW? */
fftw_complex *bgy3d_fft_axpby (DA da, fftw_complex *restrict y,
                               double alpha, double beta,
                               const fftw_complex *x)
{
    int z[3], n[3], N;

    /* Get local portion of the grid */
    DAGetCorners(da, &z[0], &z[1], &z[2], &n[0], &n[1], &n[2]);

    N = n[0] * n[1] * n[2];

    if (beta != 0.0)            /* ignore junk in y */
        for (int i = 0; i < N; i++) {
            y[i].re = alpha * x[i].re + beta * y[i].re;
            y[i].im = alpha * x[i].im + beta * y[i].im;
        }
    else
        for (int i = 0; i < N; i++) {
            y[i].re = alpha * x[i].re;
            y[i].im = alpha * x[i].im;
        }

    return y;
}

/* FIXME: constant is double so far: */
fftw_complex *bgy3d_fft_set (DA da, fftw_complex *y, double alpha)
{
    int z[3], n[3], N;

    /* Get local portion of the grid */
    DAGetCorners(da, &z[0], &z[1], &z[2], &n[0], &n[1], &n[2]);

    N = n[0] * n[1] * n[2];

    for (int i = 0; i < N; i++) {
        y[i].re = alpha;
        y[i].im = 0.0;
    }

    return y;
}


FFT_DATA *ComputeFFTfromVec(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft)

{
    int x[3], n[3];

    /* Get local portion of the grid */
    DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

    if(g_fft==NULL)
        g_fft = (FFT_DATA*) calloc(n[0] * n[1] * n[2], sizeof(*g_fft));

    /* Real Vec into complex array: */
    unpack (da, g, g_fft);

    /* forward fft */
    fft_3d(g_fft, g_fft, 1, fft_plan);

    return g_fft;
}


void ComputeVecfromFFT(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft)

{
    assert (g_fft != NULL);

    /* backward fft */
    fft_3d(g_fft, g_fft, -1, fft_plan);

    /* Pack a  complex vector with (hopefully)  vanishing imaginary part
       into a real Vec: */
    pack (da, g, g_fft);
}
