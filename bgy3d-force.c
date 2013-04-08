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


/* Long range pair potential Vec uc_fft is intent(out) here: */
static void coulomb_long_fft (const State *BHD, real G, Vec uc_fft)
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real L = PD->interval[1] - PD->interval[0];

  /* Get local portion of the k-grid */
  int x[3], n[3], i[3];
  DAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***uc_fft_;
  DAVecGetArray (BHD->dc, uc_fft, &uc_fft_);

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
            Potential,  a  complex number  with  zero imaginary  part.
            FIXME:  we  take a  square  root  here,  but in  the  most
            interesting case the Coulomb long is ~k^2 anyway:
          */
          uc_fft_[i[2]][i[1]][i[0]] =
            coulomb_long_fourier ((2 * M_PI / L) * sqrt (k2), 1.0, G);
        }
  DAVecRestoreArray (BHD->dc, uc_fft, &uc_fft_);

  /*
    Translate  the Coulomb long  field uc_fft  so that  the real-space
    representation  is localized at  the grid  center like  other grid
    representations  and not  at  the grid  corner.   This amounts  to
    scaling the k-components by a phase factor:
  */
  bgy3d_vec_fft_trans (BHD->dc, N, uc_fft);
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
  DAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***uc_fft_, ***fc_fft_[3];
  DAVecGetArray (BHD->dc, uc_fft, &uc_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fc_fft[dim], &fc_fft_[dim]);

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
  DAVecRestoreArray (BHD->dc, uc_fft, &uc_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fc_fft[dim], &fc_fft_[dim]);
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


/* Computes a pair potential. See also bgy3d_force(). */
void bgy3d_pair_potential (const DA da, const ProblemData *PD,
                           const Site a, const Site b, /* by value? */
                           Vec pot)
{
  const real G = G_COULOMB_INVERSE_RANGE;

  /* Pair   interaction  parameters:  geometric   average,  arithmetic
     average, and charge product. Site coordinates are not used. */
  const real epsilon = sqrt (a.epsilon * b.epsilon);
  const real sigma = 0.5 * (a.sigma + b.sigma);
  const real q2 = a.charge * b.charge;

  real interval[2];
  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  real h[3];
  FOR_DIM
    h[dim] = PD->h[dim];

  VecSet (pot, 0.0);

  real ***pot_;
  DAVecGetArray (da, pot, &pot_);

  int n[3], x[3], i[3];
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          real r[3];
          FOR_DIM
            r[dim] = i[dim] * h[dim] + interval[0];

          const real r_s = sqrt (SQR (r[0]) + SQR (r[1]) + SQR (r[2]));

          pot_[i[2]][i[1]][i[0]] +=
            lennard_jones_coulomb_short (r_s, sigma, epsilon, G, q2);
        }
  DAVecRestoreArray (da, pot, &pot_);
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

  const ProblemData *PD = BHD->PD;
  const DA da = BHD->da;

  real h[3];
  FOR_DIM
    h[dim] = PD->h[dim];

  const real L = PD->interval[1] - PD->interval[0];
  const real off = PD->interval[0];
  const real beta = PD->beta;

  /* FIXME: only periodic[0] == {0.0, 0.0, 0.0} appears to be used: */
  real periodic[27][3] = \
    {{0, 0, 0},
     {L, 0, 0}, {-L, 0, 0}, {0, L, 0}, {0, -L, 0}, {0, 0, L}, {0, 0, -L},
     {L, L, 0}, {-L, L, 0}, {L, -L, 0}, {-L, -L, 0},
     {L, L, L}, {-L, L, L}, {L, -L, L}, {-L, -L, L},
     {L, L, -L}, {-L, L, -L}, {L, -L, -L}, {-L, -L, -L},
     {0, L, L}, {0, -L, L}, {0, L, -L}, {0, -L, -L},
     {L, 0, L}, {-L, 0, L}, {L, 0, -L}, {-L, 0, -L}};

  FOR_DIM
    {
      VecSet (f_short[dim], 0.0);
      VecSet (f_long[dim], 0.0);
    }

  if (u_ini)
    VecSet (u_ini, 0.0);

  /* FIXME: see how c2 is being assigned, not incremented! */

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
      f_long_fft[dim] = bgy3d_vec_pop (BHD->dc);

    ComputeFFTfromCoulomb (BHD, u2, f_long, u2_fft, f_long_fft, q2);

    FOR_DIM
      bgy3d_vec_push (BHD->dc, &f_long_fft[dim]);
  }

  /*
    Sort-range  potential/force is  specific  for each  pair, on  the
    other hand:
  */
  PetscScalar ***u_ini_;
  if (u_ini)
    DAVecGetArray (da, u_ini, &u_ini_);

  PetscScalar ***c2_;
  if (c2)
    DAVecGetArray (da, c2, &c2_);

  PetscScalar ***(f_short_[3]);
  FOR_DIM
    DAVecGetArray (da, f_short[dim], &f_short_[dim]);

  /* Get local portion of the grid */
  int x[3], n[3], i[3];
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        for (int cell = 0; cell < 1; cell++) /* FIXME: one of 27 unit cells */
          {
            /* set force vectors */
            real r[3];
            FOR_DIM
              r[dim] = i[dim] * h[dim] + off + periodic[cell][dim];

            const real r_s = sqrt (SQR (r[0]) + SQR (r[1]) + SQR (r[2]));

            /*
              Lennard-Jones and Coulomb short potential. Note that the
              strength of LJ and Coulomb contributions has been scaled
              by the respective factors:
            */
            if (u_ini)
              u_ini_[i[2]][i[1]][i[0]] +=                               \
                lennard_jones_coulomb_short (r_s, sigma, epsilon, G, q2);

            /* Lennard-Jones and Coulomb short forces: */
            FOR_DIM
              f_short_[dim][i[2]][i[1]][i[0]] +=                        \
              lennard_jones_coulomb_short_grad (r_s, r[dim], sigma, epsilon, G, q2);

            /*
              Deterministic  correction.   Original  version  did  not
              respect the  damp factors for  this.  FIXME: assignment,
              not increment here as a sum over cells would imply, why?
            */
            if (c2)
              c2_[i[2]][i[1]][i[0]] =
                exp (-beta * lennard_jones_repulsive (r_s, epsilon, sigma));
          }

  if (u_ini)
    DAVecRestoreArray (da, u_ini, &u_ini_);
  if (c2)
    DAVecRestoreArray (da, c2, &c2_);
  FOR_DIM
    DAVecRestoreArray (da, f_short[dim], &f_short_[dim]);
}
