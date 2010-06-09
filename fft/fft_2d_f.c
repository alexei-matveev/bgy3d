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

#include "stdio.h"
#include "mpi.h"

#include "pack_2d.h"
#include "remap_2d.h"
#include "fft_2d.h"

/* ------------------------------------------------------------------- */
/* F77 wrapper on fft */

void fft_2d_(FFT_DATA *in, FFT_DATA *out,
	     int *flag, struct fft_plan_2d **plan)

{
  fft_2d(in,out,*flag,*plan);
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on fft_create_plan */

void fft_2d_create_plan_(
       MPI_Comm *comm, int *nfast, int *nslow,
       int *in_ilo, int *in_ihi, int *in_jlo, int *in_jhi,
       int *out_ilo, int *out_ihi, int *out_jlo, int *out_jhi,
       int *scaled, int *permute, int *nbuf,
       struct fft_plan_2d **plan)

{
  int me;

/* convert F77 indices into C */

  *plan = fft_2d_create_plan(*comm,*nfast,*nslow,
			     *in_ilo-1,*in_ihi-1,*in_jlo-1,*in_jhi-1,
			     *out_ilo-1,*out_ihi-1,*out_jlo-1,*out_jhi-1,
			     *scaled, *permute, nbuf);
  if (*plan == NULL) {
    MPI_Comm_rank(*comm,&me);
    printf("ERROR: FFT 2d plan is NULL on proc %d\n",me);
  }
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on fft_destroy_plan */

void fft_2d_destroy_plan_(struct fft_plan_2d **plan)

{
  fft_2d_destroy_plan(*plan);
}
