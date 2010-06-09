/* parallel FFT functions - 1998, 1999

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs
   (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

/* User-settable FFT precision */

/* FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag) */
/* FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag) */

#define FFT_PRECISION 2

/* ------------------------------------------------------------------------- */

/* Data types for single-precision complex */

#if FFT_PRECISION == 1

#ifdef FFT_SGI
#include "fft.h"
typedef complex FFT_DATA;
#define FFT_1D cfft1d
#define FFT_1D_INIT cfft1di
#endif

#ifdef FFT_INTEL
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft1d_
#define FFT_1D_INIT cfft1d_
#endif

#ifdef FFT_DEC
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft_
#endif

#ifdef FFT_T3E
#include <complex.h>
typedef complex single FFT_DATA;
#define FFT_1D GGFFT
#define FFT_1D_INIT GGFFT
#endif

#ifdef FFT_FFTW
#include "fftw.h"
typedef FFTW_COMPLEX FFT_DATA;
#endif

#ifdef FFT_NONE
typedef struct {
  float re;
  float im;
} FFT_DATA;
#endif

#endif

/* ------------------------------------------------------------------------- */

/* Data types for double-precision complex */

#if FFT_PRECISION == 2

#ifdef FFT_SGI
#include "fft.h"
typedef zomplex FFT_DATA;
#define FFT_1D zfft1d
#define FFT_1D_INIT zfft1di
#endif

#ifdef FFT_INTEL
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft1d_
#define FFT_1D_INIT zfft1d_
#endif

#ifdef FFT_DEC
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft_
#endif

#ifdef FFT_T3E
#include <complex.h>
typedef complex double FFT_DATA;
#define FFT_1D CCFFT
#define FFT_1D_INIT CCFFT
#endif

#ifdef FFT_FFTW
#include "fftw.h"
typedef FFTW_COMPLEX FFT_DATA;
#endif

#ifdef FFT_NONE
typedef struct {
  double re;
  double im;
} FFT_DATA;
#endif

#endif

/* ------------------------------------------------------------------------- */

/* details of how to do a 2d FFT */

struct fft_plan_2d {
  struct remap_plan_2d *pre_plan;       /* remap from input -> 1st FFTs */
  struct remap_plan_2d *mid_plan;       /* remap from 1st -> 2nd FFTs */
  struct remap_plan_2d *post_plan;      /* remap from 2nd FFTs -> output */
  FFT_DATA *copy;                   /* memory for remap results (if needed) */
  FFT_DATA *scratch;                /* scratch space for remaps */
  int total1,total2;                /* # of 1st and 2nd FFTs (times length) */
  int length1,length2;              /* length of 1st and 2nd FFTs */
  int pre_target,mid_target;        /* where to put remap results */
  int scaled;                       /* whether to scale FFT results */
  int normnum;                      /* # of values to rescale */
  double norm;                      /* normalization factor for rescaling */
                                    /* system specific 1d FFT info */
#ifdef FFT_SGI
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
#endif
#ifdef FFT_INTEL
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
#endif
#ifdef FFT_T3E
  double *coeff1;
  double *coeff2;
  double *work1;
  double *work2;
#endif
#ifdef FFT_FFTW
  fftw_plan plan_fast_forward;
  fftw_plan plan_fast_backward;
  fftw_plan plan_slow_forward;
  fftw_plan plan_slow_backward;
#endif
};

/* function prototypes */

void fft_2d(FFT_DATA *, FFT_DATA *, int, struct fft_plan_2d *);
struct fft_plan_2d *fft_2d_create_plan(MPI_Comm, int, int,
  int, int, int, int, int, int, int, int,
  int, int, int *);
void fft_2d_destroy_plan(struct fft_plan_2d *);
void factor(int, int *, int *);

/* machine specifics */

#ifdef T3E_KLUDGE

#define fft_2d_ FFT_2D
#define fft_2d_create_plan_ FFT_2D_CREATE_PLAN
#define fft_2d_destroy_plan_ FFT_2D_DESTROY_PLAN

#endif
