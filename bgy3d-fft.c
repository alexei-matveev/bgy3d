/*==========================================================*/
/*  $Id: bgy3dmolecule.c,v 1.15 2007-04-23 16:55:13 jager Exp $ */
/*==========================================================*/

#include <fftw_mpi.h>
#include "petscda.h" // DA, Vec
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

/*
 * The function  has a feature,  if the output  array g_fft is  NULL a
 * fresh array is allocated with the size derived from the distributed
 * array  description.  This  size  is suffucient  to  hold the  local
 * portion of the array.
 */
fftw_complex *ComputeFFTfromVec_fftw (DA da, fftwnd_mpi_plan fft_plan, Vec g,
				fftw_complex *g_fft, fftw_complex *scratch)
{
  int index, i[3];
  int x[3], n[3];
  PetscScalar ***g_vec;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if (g_fft == NULL)
      g_fft = (fftw_complex*) malloc(n[0] * n[1] * n[2] * sizeof(*g_fft));

  DAVecGetArray(da, g, &g_vec);
  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index = 0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  g_fft[index].re = g_vec[i[2]][i[1]][i[0]];
	  g_fft[index].im = 0;
	  index++;
	}
  DAVecRestoreArray(da, g, &g_vec);
  /* forward fft */
  fftwnd_mpi( fft_plan, 1, g_fft, scratch, FFTW_NORMAL_ORDER);

  return g_fft;
}


void ComputeVecfromFFT_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
			    fftw_complex *g_fft, fftw_complex *scratch)
{
  int index, i[3];
  int x[3], n[3];
  PetscScalar ***g_vec;

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  if(g_fft==NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Error: g_fft==NULL!\n");
      exit(1);
    }

  /* backward fft */
  fftwnd_mpi( fft_plan, 1, g_fft, scratch, FFTW_NORMAL_ORDER);

  DAVecGetArray(da, g, &g_vec);
  /* loop over local portion of grid */
  /* Attention: order of indices is not variable */
  index = 0;
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  g_vec[i[2]][i[1]][i[0]] = g_fft[index].re;
	  index++;
	}
  DAVecRestoreArray(da, g, &g_vec);


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
