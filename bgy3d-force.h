/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/* Was ealrier  #defined in bgy3d.h.  Used only  here, #undef'ed below
   after last use: */
#define CUTOFF 1.0e+8


/* Computes a pair potential. See also bgy3d_force(). */
void bgy3d_pair_potential (const State *BHD,
                           const Site a, const Site b, /* by value? */
                           Vec v_short,                /* out, real */
                           Vec v_long_fft); /* out, complex */

void bgy3d_force (State *BHD,
                  const Site a, const Site b, /* struct by value? */
                  Vec f_short[3], Vec f_long[3],
                  Vec u_ini, Vec c2,
                  Vec u2, Vec u2_fft,
                  real damp, real damp_LJ);

static inline real
lennard_jones (real r, real epsilon, real sigma)
{
  if (unlikely ((sigma == 0.0 || epsilon == 0.0) && r == 0.0))
    return 0.0;                 /* 0/0 == nan */
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
  if (unlikely ((sigma == 0.0 || epsilon == 0.0) && r == 0.0))
    return 0.0;                 /* 0/0 == nan */
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

    coulomb_long_fourier (k, 1.0, G) = (4π/ε₀) exp (- k² / 4G²) / k²

  Note that the long-range Coulomb is regular in real space:

                            2           4
    erf(x) / x = 2 / √π - 2x / 3√π + O(x )

  So that for a typical value of G = 1.2 the long range potential at r
  = 0  is 331.84164 * 1.2 *  1.12837916709551 ~ 449 kcal  for two unit
  charges. By doing a finite-size FFT of the k-gitter appriximation of
  the  above  Fourier representation  you  will  likely get  something
  different.
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


/*
  FIXME: Functions marked as deprecated are actually unused. They also
  should not be used as  long-range interactions are best specified on
  the k-grid and not on the r-grid.
*/
static inline real deprecated
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

static inline real deprecated
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

static inline real deprecated
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

static inline real deprecated
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

#undef CUTOFF


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
