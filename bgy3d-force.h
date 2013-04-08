/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O.c,v 1.42 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

// #include "bgy3d.h"              /* real, EPSILON0INV, CUTOFF, SQR() */

/* Computes a pair potential. See also bgy3d_force(). */
void bgy3d_pair_potential (const DA da, const ProblemData *PD,
                           const Site a, const Site b, /* by value? */
                           Vec pot);

void bgy3d_force (State *BHD,
                  const Site a, const Site b, /* struct by value? */
                  Vec f_short[3], Vec f_long[3],
                  Vec u_ini, Vec c2,
                  Vec u2, Vec u2_fft,
                  real damp, real damp_LJ);

static inline real
lennard_jones (real r, real epsilon, real sigma)
{
  const real sr = sigma / r;
  const real sr6 = SQR (sr) * SQR (sr) * SQR (sr);

  const real re = 4 * epsilon * sr6 * (sr6 - 1);

  if (fabs (re) > epsilon * CUTOFF)
    return epsilon * CUTOFF;
  else
    return re;
}

static inline real
lennard_jones_grad (real r, real xr, real epsilon, real sigma)
{
  if (xr == 0.0)
    return 0.0;
  if (r == 0.0)
    return -epsilon * CUTOFF;

  const real sr = sigma / r;
  const real sr6 = SQR (sr) * SQR (sr) * SQR (sr);

  const real re = -24 * epsilon * sr6 / r * (2 * sr6 - 1) * xr / r;

  if (re > epsilon * CUTOFF)
    return epsilon * CUTOFF;
  else if (re < -epsilon * CUTOFF)
    return -epsilon * CUTOFF;
  else
    return re;
}

/* This  one  is only  used  for  some  kind of  pre-conditioning,  or
   regularization inside the solvent excluded volume: */
static inline real
lennard_jones_repulsive (real r, real epsilon, real sigma)
{
  const real sr = sigma / r;
  const real sr6 = SQR (sr) * SQR (sr) * SQR (sr);

  const real re = 4 * epsilon * sr6 * sr6;

  if (fabs (re) > epsilon * CUTOFF)
    return epsilon * CUTOFF;
  else
    return re;
}

/* NOTE: so far  in all cases the returned  result contains the factor
   q2. */
static inline real
coulomb_short (real r, real q2, real G)
{
  if (r == 0.0)
    return EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
  else
    {
      real re = EPSILON0INV * q2 * erfc (G * r) / r;

      /* Check for large values remember: exp(-re) will be computed: */
      if (fabs (re) > fabs (EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
        return EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
      else
        return re;
    }
}

static inline real
coulomb_short_grad (real r, real rx, real q2, real G)
{
  if (rx == 0)
    return 0.0;

  if (r == 0.0)
    return -EPSILON0INV * q2 * (CUTOFF*1.0e-5);
  else
    {
      real re = - EPSILON0INV * q2 * (erfc (G * r) +
                                      2. * G / sqrt (M_PI) * r * exp(-G * G * r * r)) * rx / pow (r, 3.0);

      if (fabs (re) > fabs (EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
        return -EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
      else
        return re;
    }
}


/*
  In the most general case

  coulomb_long (r, 1.0, G) = (1/ε₀) erf (G r) / r

  An  extra   argument  q2  is   the  overall  scaling   factor.   The
  corresponding Fourier transform is

  coulomb_long_Fourier (k, 1.0, G) = (4π/ε₀) exp (- k² / 4G²)
*/
static inline real
coulomb_long_fourier (real k, real q2, real G)
{
  const real k2 = SQR (k);

  if (unlikely (k2 == 0.0))
    return 0.0;         /* 1/k2 is undefined */
  else
    return (4 * M_PI * EPSILON0INV) * q2 * exp (-k2 / (4 * SQR (G))) / k2;
}


static inline real
coulomb_long (real r, real q2, real G)
{
   if (r == 0.0)
     return EPSILON0INV * q2 * G * 2.0 / sqrt(M_PI);
   else
     {
       real re = EPSILON0INV * q2 * erf (G * r) / r;

       /* Check for large values, remember: exp(-re) will be computed */
       if (fabs (re) > fabs (EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
         return EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
       else
         return re;
     }
}

static inline real
coulomb_long_grad (real r, real rx, real q2, real G)
{
  if (r == 0.0)
    return 0.0;
  else
    {
      real re = - EPSILON0INV * q2 * (erf (G * r)
                                      - 2. * G / sqrt (M_PI) * r * exp(-G * G * r * r)) * rx / pow(r,3.0);

      if (fabs (re) > fabs (EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
        return -EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
      else
        return re;
    }
}

static inline real
coulomb (real r, real q2)
{
   if (r == 0.0)
     return EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
   else
     {
       real re = EPSILON0INV * q2 /r;

       if (fabs (re) > fabs (EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
         return EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
       else
         return re;
     }
}

static inline real
coulomb_grad (real r, real rx, real q2)
{
  if (rx == 0)
    return 0;

  if (r == 0)
    return -EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
  else
    {
      real re = - EPSILON0INV * q2 * rx / pow (r, 3.0);

      if (fabs (re) > fabs (EPSILON0INV * q2 * (CUTOFF * 1.0e-5)))
        return -EPSILON0INV * q2 * (CUTOFF * 1.0e-5);
      else
        return re;
    }
}


/*
  Short-range part of  the interaction of two charged  LJ site (sigma,
  epsilon, charge).  Note that for  pair interaction one  will usually
  pass the  arithmetic average of  two sigmas, a geometric  average of
  two epsilons and the charge product as arguments
 */
static inline real
lennard_jones_coulomb_short (real r,
                             real sigma, real epsilon,
                             real G, real charge)
{
  /*
    1.   Lennard-Jones  and  2.    Coulomb,  short  range  part.   For
    historical reasons  the overall scaling factors,  e.g. the product
    of  solvent-  and  solute  site  charges  q,  is  handled  by  the
    elementary functions:
  */
  return
    lennard_jones (r, epsilon, sigma) +
    coulomb_short (r, charge, G);
}


static inline real
lennard_jones_coulomb_short_grad (real r, real rx,
                                  real sigma, real epsilon,
                                  real G, real charge)
{
  return
    lennard_jones_grad (r, rx, epsilon, sigma) +
    coulomb_short_grad (r, rx, charge, G);
}
