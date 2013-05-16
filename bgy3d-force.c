/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O.c,v 1.42 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* G_COULOMB_INVERSE_RANGE, needs Site */
#include "bgy3d-vec.h"          /* bgy3d_vec_map*() */
#include <complex.h>            /* after fftw.h */
#include "bgy3d-force.h"        /* Coulomb_short() etc. */

/* FIXME: bgy3d-solvents.h pollutes the namespace: */
#undef sH
#undef eH
#undef qH
#undef sO
#undef eO
#undef qO


/* Tabulate v = f(r) with origin at the grid center:  */
static void vec_rmap (const State *BHD, real (*f)(real r), Vec v)
{
  const real off = BHD->PD->interval[0];

  real h[3];
  FOR_DIM
    h[dim] = BHD->PD->h[dim];

  real ***v_;
  DMDAVecGetArray (BHD->da, v, &v_);

  int n[3], x[3], i[3];
  DMDAGetCorners (BHD->da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          real r[3];
          FOR_DIM
            r[dim] = i[dim] * h[dim] + off;

          const real r_s = sqrt (SQR (r[0]) + SQR (r[1]) + SQR (r[2]));

          v_[i[2]][i[1]][i[0]] = f (r_s);
        }
  DMDAVecRestoreArray (BHD->da, v, &v_);
}


/* Tabulate v_fft = f(k) with origin at the grid corner: */
static void vec_kmap (const State *BHD, complex (*f)(real k), Vec v_fft)
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real L = PD->interval[1] - PD->interval[0];

  /* Get local portion of the k-grid */
  int x[3], n[3], i[3];
  DMDAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***v_fft_;
  DMDAVecGetArray (BHD->dc, v_fft, &v_fft_);

   /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);

          /*
            FIXME:  we  take a  square  root  here,  but in  the  most
            interesting case  the Coulomb depends on  k^2 anyway. This
            potential  is also  a complex  number with  zero imaginary
            part.
          */
          v_fft_[i[2]][i[1]][i[0]] = f ((2 * M_PI / L) * sqrt (k2));
        }
  DMDAVecRestoreArray (BHD->dc, v_fft, &v_fft_);
}


/* Long  range  pair potential  Vec  uc_fft  is intent(out),  complex,
   centered Vec: */
static void coulomb_long_fft (const State *BHD, real G, Vec uc_fft)
{
  /*
    Tabulate  spherically  symmetric   function  around  grid  corner.
    Potential is a complex number with zero imaginary part.
  */
  complex pure f (real k)
  {
    return coulomb_long_fourier (k, 1.0, G);
  }
  vec_kmap (BHD, f, uc_fft);

  /*
    Translate  the Coulomb long  field uc_fft  so that  the real-space
    representation  is localized at  the grid  center like  other grid
    representations  and not  at  the grid  corner.   This amounts  to
    scaling the k-components by a phase factor:
  */
  bgy3d_vec_fft_trans (BHD->dc, BHD->PD->N, uc_fft);
}


/* Computes FFT of the gradients  from FFT of a scalar.  FIXME: rather
   these are forces! */
static void grad_fft (const State *BHD, Vec uc_fft, Vec fc_fft[3])
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real L = PD->interval[1] - PD->interval[0];
  const real fac = 2.0 * M_PI / L;

  /* Get local portion of the k-grid */
  int x[3], n[3], i[3];
  DMDAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***uc_fft_, ***fc_fft_[3];
  DMDAVecGetArray (BHD->dc, uc_fft, &uc_fft_);
  FOR_DIM
    DMDAVecGetArray (BHD->dc, fc_fft[dim], &fc_fft_[dim]);

   /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          /* Force  is purely imaginary  if Vec  uc_fft happens  to be
             real: */
          FOR_DIM
            fc_fft_[dim][i[2]][i[1]][i[0]] = ic[dim] * ((fac * I) * uc_fft_[i[2]][i[1]][i[0]]);
        }
  DMDAVecRestoreArray (BHD->dc, uc_fft, &uc_fft_);
  FOR_DIM
    DMDAVecRestoreArray (BHD->dc, fc_fft[dim], &fc_fft_[dim]);
}


/*
  Long range  pair potential Vec uc  is intent(out) here,  same as its
  FFT  transform uc_fft.  Vec  fc[] is  filled with  the corresponding
  force derived by means of FFT from the potential uc.

  No side effects.
*/
static void ComputeFFTfromCoulomb (State *BHD,
                                   Vec uc, Vec fc[3], /* intent(out) */
                                   Vec uc_fft,    /* complex, intent(out) */
                                   Vec fc_fft[3], /* complex, intent(out) */
                                   real factor)
{
  const real G = G_COULOMB_INVERSE_RANGE;

  /* Potential of a unit charge located at the grid center: */
  coulomb_long_fft (BHD, G, uc_fft);

  VecScale (uc_fft, factor);

  /* Corresponding forces: */
  grad_fft (BHD, uc_fft, fc_fft);

  const ProblemData *PD = BHD->PD;
  const real L = PD->interval[1] - PD->interval[0];

  /* FFT^-1 potential ... */
  MatMultTranspose (BHD->fft_mat, uc_fft, uc);
  VecScale (uc, 1./L/L/L);

  /* FFT^-1 of the corresponding force: */
  FOR_DIM
    {
      MatMultTranspose (BHD->fft_mat, fc_fft[dim], fc[dim]);
      VecScale (fc[dim], 1./L/L/L);
    }
}


/* Computes a pair potential. Centered. See also bgy3d_force(). */
void bgy3d_pair_potential (const State *BHD,
                           const Site a, const Site b, /* by value? */
                           Vec v_short,    /* out, real, center */
                           Vec v_long_fft) /* out, complex, center */
{
  const real G = G_COULOMB_INVERSE_RANGE;

  /* Pair   interaction  parameters:  geometric   average,  arithmetic
     average, and charge product. Site coordinates are not used. */
  const real epsilon = sqrt (a.epsilon * b.epsilon);
  const real sigma = 0.5 * (a.sigma + b.sigma);
  const real q2 = a.charge * b.charge;

  /* Tabulate spherically symmetric function around grid center: */
  real pure f (real r)
  {
    return lennard_jones_coulomb_short (r, sigma, epsilon, G, q2);
  }
  vec_rmap (BHD, f, v_short);

  /* Long-range part of the potential is best represented by FFT: */
  coulomb_long_fft (BHD, G, v_long_fft);
  VecScale (v_long_fft, q2);
}


/* Precompute forces and more for a pair: */
void bgy3d_force (State *BHD,
                  const Site a, const Site b, /* struct by value? */
                  Vec f_short[3], Vec f_long[3],
                  Vec u_ini, Vec c2,
                  Vec u2, Vec u2_fft,
                  real damp, real damp_LJ)
{
  const real G = G_COULOMB_INVERSE_RANGE;
  /*
    Pair   interaction  parameters:   arithmetic   average,  geometric
    average, and  charge product.   Parameters that define  the energy
    scale are scaled by damping factors.
  */
  const real sigma = 0.5 * (a.sigma + b.sigma);
  const real epsilon = damp_LJ * sqrt (a.epsilon * b.epsilon);
  const real q2 = damp * (a.charge * b.charge);

  const real beta = BHD->PD->beta;

  /*
    Compute long-range  Coulomb potential and  corresponding forces by
    FFT.

    Here Vec u2  and a complex array u2_fft[]  both are intent(out) in
    the next call.  The Vec f_long[], intent(out), is  filled with the
    corresponding force.   Performs 4 FFTs.  Again note that  the only
    difference for all u2[i][j] and their FFT transform is the overall
    scaling factor  q[i] * q[j].  FIXME: why  keeping O(m^2) versions,
    with m being number of solvent sites, of almost the same field and
    repeating unnecessary FFTs?

    NOTE: this code fills temp  Vec f_long_fft[3] with FFT of the long
    range Coulomb forces and discards that data.
  */
  {
    local Vec f_long_fft[3];
    FOR_DIM
      f_long_fft[dim] = vec_pop (BHD->dc);

    ComputeFFTfromCoulomb (BHD, u2, f_long, u2_fft, f_long_fft, q2);

    FOR_DIM
      vec_push (BHD->dc, &f_long_fft[dim]);
  }

  /*
    Sort-range  potential/force is  specific  for each  pair, on  the
    other hand.

    Lennard-Jones and Coulomb short  potential. Note that the strength
    of LJ and Coulomb contributions  has been scaled by the respective
    factors:
  */
  if (u_ini)
    {
      real pure f (real r)
      {
        return lennard_jones_coulomb_short (r, sigma, epsilon, G, q2);
      }
      vec_rmap (BHD, f, u_ini);
    }

  /*
    Deterministic  correction.  Original version  did not  respect the
    damp factors for this.
  */
  if (c2)
    {
      real pure f (real r)
      {
        /* FIXME: argument order? */
        return exp (-beta * lennard_jones_repulsive (r, epsilon, sigma));
      }
      vec_rmap (BHD, f, c2);
    }


  /* Now tabulate the short-range forces: */
  real h[3];
  FOR_DIM
    h[dim] = BHD->PD->h[dim];

  const real off = BHD->PD->interval[0];

  PetscScalar ***(f_short_[3]);
  FOR_DIM
    DMDAVecGetArray (BHD->da, f_short[dim], &f_short_[dim]);

  /* Get local portion of the grid */
  int x[3], n[3], i[3];
  DMDAGetCorners (BHD->da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
          {
            /* set force vectors */
            real r[3];
            FOR_DIM
              r[dim] = i[dim] * h[dim] + off;

            const real r_s = sqrt (SQR (r[0]) + SQR (r[1]) + SQR (r[2]));

            /* Lennard-Jones and Coulomb short forces: */
            FOR_DIM
              f_short_[dim][i[2]][i[1]][i[0]] =
                lennard_jones_coulomb_short_grad (r_s, r[dim], sigma, epsilon, G, q2);
          }

  FOR_DIM
    DMDAVecRestoreArray (BHD->da, f_short[dim], &f_short_[dim]);
}
