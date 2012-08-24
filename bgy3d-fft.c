/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.15 2007-04-23 16:55:13 jager Exp $ */
/*==========================================================*/

#include <assert.h>
#include <fftw_mpi.h>
#include "fft_3d.h"             /* FFT_DATA */
#include "petscda.h"            /* DA, Vec */
#include "bgy3d-fft.h"

fftw_complex *bgy3d_fft_malloc (DA da)
{
    int x[3], n[3];

    /* Get local portion of the grid */
    DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

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
    DAGetCorners(da, &(i0), &(j0), &(k0), &(ni), &(nj), &(nk));

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
    DAGetCorners(da, &(i0), &(j0), &(k0), &(ni), &(nj), &(nk));

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
    DAGetCorners(da, &(z[0]), &(z[1]), &(z[2]), &(n[0]), &(n[1]), &(n[2]));

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
    DAGetCorners(da, &(z[0]), &(z[1]), &(z[2]), &(n[0]), &(n[1]), &(n[2]));

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
    DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

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
