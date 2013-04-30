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
#include "bgy3d-vec.h"          /* bgy3d_vec_destroy() */
#include "bgy3d-poisson.h"
#include <complex.h>            /* after fftw.h */

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
  Vec work = bgy3d_vec_create (BHD->dc);

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
    int x[3], n[3], i[3];
    DMDAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

    complex ***work_;
    DMDAVecGetArray (BHD->dc, work, &work_);

    for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
      for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
        for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
          {
            /* FIXME: what if we  change the complex vectors to remove
               the redundnacy? */
            int ic[3];

            /* Take negative frequencies for i > N/2: */
            FOR_DIM
              ic[dim] = KFREQ (i[dim], N[dim]);

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
    DMDAVecRestoreArray (BHD->dc, work, &work_);
  }

  /* u(x, y, z) := IFFT(u(kx, ky, kz)) */
  MatMultTranspose (BHD->fft_mat, work, uc);

  bgy3d_vec_destroy (&work);
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
  Vec work = bgy3d_vec_create (BHD->dc);

  /* Get FFT of  rho: rho(i, j, k) -> fft_rho(kx,  ky, kz) placed into
     complex work: */
  MatMult (BHD->fft_mat, rho, work);

  /* With q = -4π/ε₀ you would get the potential: */
  const real scale = - q / NNN / 4;

  /* Loop over local portion of the k-grid */
  {
    int i0, j0, k0, ni, nj, nk;
    DMDAGetCorners (BHD->dc, &i0, &j0, &k0, &ni, &nj, &nk);

    complex ***work_;
    DMDAVecGetArray (BHD->dc, work, &work_);

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
    DMDAVecRestoreArray (BHD->dc, work, &work_);
  }

  /* u(x, y, z) := IFFT(u(kx, ky, kz)) */
  MatMultTranspose (BHD->fft_mat, work, uc);

  bgy3d_vec_destroy (&work);
}
#endif
