/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li

  The  funciton  bgy3d_solute_field()  initiates  the  computation  of
  (short  range)  site-specifc  force  field and  long  range  Coulomb
  potential (common  for all sites)  for arbitrary number  of solvent-
  and solute sites.
*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-poisson.h"      /* bgy3d_poisson() */
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* G_COULOMB_INVERSE_RANGE */
#include "bgy3d-force.h"        /* bgy3d_coulomb_long_fft() */
#include "bgy3d-vec.h"          /* bgy3d_vec_read() */


/* Singular as 4/r¹² */
static inline pure real
lj0 (real r)
{
  const real p6 = 1 / pow (r, 6);
  return 4 * p6 * (p6 - 1);
}


/* Singular as -48/r¹³ */
static inline pure real
lj1 (real r)
{
  const real p6 = 1 / pow (r, 6);
  return 4 * p6 * (6 - 12 * p6) / r;
}


#define RMIN 0.02               /* 2% of sigma */
/*
  Regularized at R,  the same as lj0(r) for r >=  R.  For lower values
  of the distance  r we reperesent potential by  a cap, an upside-down
  parabola  a + b  * (R²  - r²).   To glue  the two  cases we  need to
  satisfy these at r = R:

    a + b(r² - R²) = f(r)
    2br = f'(r)
*/
static inline pure real
ljcap0 (real r)
{
  const real R = RMIN;
  if (likely (r >= R))
    return lj0 (r);
  else
    {
      // printf ("ljcap0: r = %f\n", r);
      const real b = lj1 (R) / (2 * R);
      const real a = lj0 (R);
      return a + b * (r - R) * (r + R);
    }
}


/* Regularized at R */
static inline pure real
ljcap1 (real r)
{
  const real R = RMIN;
  if (likely (r >= R))
    return lj1 (r);
  else
    {
      // printf ("ljcap1: r = %f\n", r);
      const real b = lj1 (R) / (2 * R);
      return  2 * b * r;
    }
}
#undef RMIN


/* Singlular as 1/r */
static inline pure real
cs0 (real r)
{
  return erfc (r) / r;
}


/* Singlular as -1/r²: */
static inline pure real
cs1 (real r)
{
  const real r2 = pow (r, 2);
  return - (erfc (r) + (2 / sqrt (M_PI)) * r * exp (-r2)) / r2;
}


#define RMIN 0.02               /* 2% of G^-1 */
/*
  Regularized at R,  the same as cs0(r) for r >=  R.  For lower values
  of the distance  r we reperesent potential by  a cap, an upside-down
  parabola  a + b  * (R²  - r²).   To glue  the two  cases we  need to
  satisfy these at r = R:

    a + b(r² - R²) = f(r)
    2br = f'(r)
*/
static inline pure real
cscap0 (real r)
{
  const real R = RMIN;
  if (likely (r >= R))
    return cs0 (r);
  else
    {
      // printf ("cscap0: r = %f\n", r);
      const real b = cs1 (R) / (2 * R);
      const real a = cs0 (R);
      return a + b * (r - R) * (r + R);
    }
}
#undef RMIN


/*
  Smooth long-range erf(r)/r, finite:

    taylor ((√π/2) * erf(r)/r, r, 0, 8) =

            2    4    6    8
           r    r    r    r
       1 - -- + -- - -- + --- + . . .
           3    10   42   216

*/
static inline pure real
cl0 (real r)
{
  if (likely (r >= 0.02))       /* 0.02^8/216 ~ 1.2e-16 */
    return erf (r) / r;
  else
    {
      // printf ("cl0: r = %f\n", r);
      const real r2 = pow (r, 2);
      return (2 / sqrt (M_PI)) *
        (1. + r2 * (-1./3. + r2 * (1./10. + r2 * (-1./42.))));
    }
}


/*
  Derivative of erf(r)/r is also finite and smooth, even if it may not
  look so:

       -r²
    2 e         erf(r)
   --------  -  ------
    √π r          r²


  Taylor ((-3√π/4r) * derivative (erf(r)/r, r), r, 0, 8) =
               2      4    6    8
            3 r    3 r    r    r
        1 - ---- + ---- - -- + -- + . . .
             5      14    18   88

*/
static inline pure real
cl1 (real r)
{
  const real r2 = pow (r, 2);
  if (likely (r >= 0.02))       /* 0.02^8/88 ~ 2.9e-16 */
    return ((2 / sqrt (M_PI)) * r * exp (-r2) - erf (r)) / r2;
  else
    {
      assert (false);
      return - (4 / (3 * sqrt (M_PI))) * r *
        (1. + r2 * (-3./5. + r2 * (3./14. + r2 * (-1./18.))));
    }
}



/*
  Calculate a  real field  "f" for the  solvent site  characterized by
  (sigma, epsilon,  charge) in the presence  of the solute  S at every
  point (x, y, z) of the local grid.

  This  function expects  a callback  f() obeying  specific interface.
  The function f(x, y,  z, eps, sig, chg, S) can be  ljc() or rho() as
  two examples.

  Site  A  is  intent(in) but  has  to  be  passed by  value  (copying
  involved) at  the moment.  We  are modifying coordinates of  a local
  copy  and pass  a reference  to that  copy further  to  the callback
  function *f.  The coordinates of  a site as passed to this functions
  are ignored.

  Vector "v" is the intent(out) argument.
 */
static void
field (const State *BHD, Site A, int n, const Site S[n],
       real (*f)(const Site *A, int n, const Site S[n]),
       Vec v)
{
  /*
    Compute the field  f at (x, y,  z) <-> (i, j, k)  e.g.  by summing
    (LJ) contributions from all solute sites at that grid point:
  */
  real f3 (const real r[3])
  {
    FOR_DIM
      A.x[dim] = r[dim];        /* Modifying the input! */
    return f (&A, n, S);
  }

  vec_rmap3 (BHD, f3, v);
}


/* Callback here is a function passed from QM code: */
static void
grid_map (DA da, const ProblemData *PD,
          void (*f)(int m, const real x[m][3], real fx[m]),
          Vec v)
{
  int i0, j0, k0;
  int ni, nj, nk;

  /* Get local portion of the grid */
  DMDAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* MEMORY: huge arrays here: */
  int m = ni * nj * nk;

  /*
    Use  dynamically  memory  allocation  here, otherwise  stack  will
    overflow.  FIXME:  do the work  in chunks if  memory/efficiency is
    low (see bgy3d-potential.c for example). This allocation is as big
    as four real Vecs:
  */
  real (*x)[3], *fx;
  fx = malloc (m * sizeof (real));
  x = malloc (m * 3 * sizeof (real)); /* contiguous memory */

  /* Get coordinates of the local grid portion: */
  {
    int ijk = 0;
    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          {
            /* Coordinates (x, y, z) <-> (i, j, k): */
            x[ijk][0] = i * PD->h[0] - PD->L[0] / 2;
            x[ijk][1] = j * PD->h[1] - PD->L[1] / 2;
            x[ijk][2] = k * PD->h[2] - PD->L[2] / 2;
            ijk++;
          }
    assert (ijk == m);
  }

  /* Let the function f() compute  the field/density at our portion of
     the grid. FIXME: cast here is to silence the const-warning: */
  f (m, (const real (*)[3]) x, fx);

  /* Copy contents of fx[] to the output vector: */
  {
    PetscScalar ***v_;
    DMDAVecGetArray (da, v, &v_);
    int ijk = 0;
    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          {
            v_[k][j][i] = fx[ijk];
            ijk++;
          }
    assert (ijk == m);
    DMDAVecRestoreArray (da, v, &v_);
  }

  /* Remember to free! */
  free (fx);
  free (x);
}


/*
  Interaction of a charged LJ  site (sigma, epsilon, charge) at (x, y,
  z) with the solute S.

  This function obeys the callback interface assumed in field(). It is
  supposed  to get  (1) parameters  of the  solvent site  such  as its
  location and  force field parameters,  and (2) a description  of the
  solute and return a real number such as an interaction energy or the
  charge density:
*/
static real
ljc (const Site *A, int n, const Site S[n])
{
  const real G = G_COULOMB_INVERSE_RANGE;

  /* Sum force field contribution from all solute sites: */
  real field = 0.0;

  for (int i = 0; i < n; i++)
    {
      const Site *B = &S[i];    /* shorter alias */

      /* Interaction parameters for a pair of LJ sites: */
      real eab = sqrt (A->epsilon * B->epsilon);
      real sab = 0.5 * (A->sigma + B->sigma);
      real qab = A->charge * B->charge * EPSILON0INV;

      /* Distance from a grid point to this site: */
      real rab = sqrt (SQR(A->x[0] - B->x[0]) +
                       SQR(A->x[1] - B->x[1]) +
                       SQR(A->x[2] - B->x[2]));

      /* Lennard-Jones + Coulomb, short range part: */
      field += eab * ljcap0 (rab / sab) + qab * G * cscap0 (G * rab);
    }

  return field;
}


/*
  Gaussian functions  represent each  solute site charge,  compute the
  total charge density at every point x[i].

  Each gaussian is evaluated as:

    ρ(r) = q * [G / √π]³ * exp[-G² * (r - x₀)²]

  NOTE: one  might want to make  the gaussian width a  property of the
  (solute) site, but consider  the implications first.  At some places
  we assume that the  long-range Coulomb potential acting on (solvent)
  sites differs only by a factor  (charge of the solvent site) to save
  space. This  would not  be the case  if solvent sites  had different
  widths.  One  might  of   course  treat  solute  and  solvent  sites
  differently, but why?
*/
static void
cores (const State *BHD,
       int n, const real q[n], real r[n][3], real G, /* in */
       Vec rho)                 /* out, real, center */
{
  real f3 (const real x[3])
  {
    /* Sum Gaussian contributions from all (solute) sites: */
    real sum = 0.0;
    for (int i = 0; i < n; i++)
      {
        /* Square of the distance from a grid point to this site: */
        real r2 = SQR (x[0] - r[i][0]) +
                  SQR (x[1] - r[i][1]) +
                  SQR (x[2] - r[i][2]);

        /* Gaussian  distribution, note  that G  is not  a  width, but
           rather an inverse of it: */
        sum += q[i] * exp(- SQR (G) * r2);
      }
    return pow (G / sqrt (M_PI), 3.0) * sum;
  }
  vec_rmap3 (BHD, f3, rho);
}


/*
  Computes the weighted "form factor"

    f(k) = Σ  q  exp (-I kx )
            i  i           i
 */
static inline complex
form_factor (const real k[3], int m, const real q[m], real x[m][3])
{
  /*
    This form factor is regular for  small k with a limit equal to the
    total charge of the system and the higher order given by multipole
    expansion:
  */
  complex sum = 0.0;
  for (int i = 0; i < m; i++)
    sum += q[i] * cexp (-I * (k[0] * x[i][0] +
                              k[1] * x[i][1] +
                              k[2] * x[i][2]));
  return sum;
}


/*
  Updates  v_fft with  its  convolution with  the "instantaneous  core
  density" of solutes. Technically, scales the k-component v(k) by the
  form   factor   derived  from   the   weights   q  and   coordinates
  x. Intuitively, replaces v(r) of a single unit charge at origin by a
  weighted superposition of such potentials  from every core in x. You
  choose which definition fits the most.
 */
void
bgy3d_solute_form (const State *BHD,
                   int m, const real q[m], real x[m][3],
                   Vec v_fft)   /* inout */
{
  local Vec f_fft = vec_duplicate (v_fft);

  /* Tabulate form factor: */
  complex pure f3 (const real k[3])
  {
    return form_factor (k, m, q, x);
  }
  vec_kmap3 (BHD, f3, f_fft);

  /* Pointwise-product in k-space: */
  complex pure mul (complex a, complex b)
  {
    return a * b;
  }
  vec_fft_map2 (v_fft, mul, f_fft, v_fft); /* aliasing! */

  vec_destroy (&f_fft);
}

/*
  In  the   spherically  symmetric  single  center   case  the  smooth
  assymptotic potential is chosen as

    u(r) = (1/ε₀) erf (G r) / r

  The corresponding Fourier transform is

    u(k) = (4π/ε₀) exp (- k² / 4G²) / k²

  Note that the long-range Coulomb is regular in real space:

    erf(x) / x = 2 / √π - 2x² / 3√π + O(x⁴)

  So that for a typical value of G = 1.2 the long range potential at r
  = 0  is 331.84164 * 1.2 *  1.12837916709551 ~ 449 kcal  for two unit
  charges.  By doing a finite  size FFT of the k-grid appriximation of
  the  above  Fourier representation  you  will  likely get  something
  different.

  For many centers we effectively  "fuse" several logical steps in one
  procedure:

  1) Tabulate electric form factor  ρ(k) for a superposition of δ-like
  functions  at  x[i] with  weights  q[i].   Centered  at the  corner.
  Normalization is such that

    ρ(k=0) = Σ q[i]
              i

  The phase factor accounts for the translation.  If this ever becomes
  performance   critical  consider   using  the   properties   of  the
  exponential and maybe constant grid spacing.

  2) Obtain ρ(k) for a  superposition of gaussians by convolution with
  a unit gaussian proportional to

    exp(-k² / 4G²)

  3) "Solve" a Poisson equation by scaling with (4π/ε₀)/k².

  4) Translate  the whole  to the grid  center by flipping  the phase.
  Translation    and   convolutions   (both    effectively   pointwise
  multiplications in k-space) are commutative, by the way.

  We re-use the single core version from bgy3d-force.c here:
*/
static void
coulomb_long_fft (const State *BHD,
                  int n, const real q[n], real x[n][3], real G, /* in */
                  Vec uc_fft)   /* out, complex, center */
{
  /*
    This is a long-range Coulomb shape around the grid center.  FIXME:
    Cannot divide by zero here  to compute something that is inversely
    proportional  to  k², but  this  is  the  value that  defines  the
    "average" potential in the cell:
  */
  bgy3d_coulomb_long_fft (BHD, G, uc_fft);

  /*
    Convolution with the charge-weighted  form factor. The form factor
    at k=0 is zero for neutral  species, so that in this case a finite
    shift  in uc(k=0)  would not  matter.  FIXME:  1/k² is  not finite
    though, hopefully we get a meaningful dipole field here:
  */
  bgy3d_solute_form (BHD, n, q, x, uc_fft);
}


static void
coulomb_long (const State *BHD,
              int n, const real q[n], real x[n][3], real G, /* in */
              Vec uc)           /* out, real, center */
{
  real f3 (const real y[3])
  {
    /* Sum contributions from all (solute) sites: */
    real sum = 0.0;
    for (int i = 0; i < n; i++)
      {
        /* Distance from a grid point to this site: */
        real r = sqrt (SQR (y[0] - x[i][0]) +
                       SQR (y[1] - x[i][1]) +
                       SQR (y[2] - x[i][2]));

        /* cl0(r) == erf(r)/r */
        sum += q[i] * G * cl0 (G * r);
      }
    return EPSILON0INV * sum;
  }
  vec_rmap3 (BHD, f3, uc);
}


/*
  Create initial solute field.

  See Eqs.  (5.106)  and (5.08) in the Jager  thesis.  Fill us[i] with
  short-range potential acting on site i,

    u  = v   + v
     S    LJ    CS

  and uc_fft with k-space representation of long range (or rather
  asymptotic) electrostatic potential

            -1
    v  = FFT  (uc_fft)
     CL

  not scaled by the charges of the solvent sites.

  Returns both, the asymptotic  electrostatic potential Vec uc_fft and
  the corresponding  (diffuse part of the) charge  density Vec uc_rho.
  If both, uc_fft and uc_rho, are NULL the asymptotic Coulomb field is
  not computed. FIXME: Otherwise both must be non-NULL.

  The  function pointer density()  passed to  bgy3d_solute_field(), if
  not NULL,  is used to  compute the density of  (negatively charged!)
  electrons at arbitrary point in space.
*/
void
bgy3d_solute_field (const State *BHD,
                    int m, const Site solvent[m], /* in */
                    int n, const Site solute[n],  /* in */
                    Vec us[m],                    /* out, optional */
                    Vec uc_fft,                   /* out, optional */
                    Vec uc_rho,                   /* out, optional */
                    Vec uc,                       /* out, optional */
                    void (*density)(int k, const real x[k][3], real rho[k]))
{
  const real G = G_COULOMB_INVERSE_RANGE;

  PRINTF ("Computing solute data\n");


  /*
    Calculate FF potential for all solvent sites.

    Scaling  the  epsilon of  all  sites by  factor  X  scales the  LJ
    site-site interaction with by the factor X.  Use this fact instead
    of  handling  damp  factors  separately.  Similarly,  scaling  the
    solvent site  charge by a  factor scales the  Coulomb interaction.
    At least in  the two test examples the  two factors are identical.
    Initial  version of the  code scaled  LJ- and  short-range Coulomb
    interaction using  two different factors.

    If want to scale all the interactions uniformly, you may choose to
    scale the output of this funciton.
  */

  /*
    Fill the force field interaction  of solvent sites with the solute
    of general shape into the respective arrays.

    We supply ljc() as a callback function that is supposed to compute
    the interaction of a charged LJ solvent site with the solute.
   */
#ifndef QM
  const real scale_coul_short = 1.0; /* regular case */
#else
  const real scale_coul_short = 0.0;
#endif
  /*
    In a special  case (ifdef QM) the short  range Coulomb interaction
    of the  solvent site with  the solute was deliberately  omitted by
    scaling charge of the solvent site down to zero.  This effectively
    makes  a  point  charge  (Coulomb  short  +  Coulomb  long)  to  a
    distributed Gaussian (Coulomb long only).
  */
  if (us)    /* Not quite sure if passing us = NULL is legal though */
    for (int i = 0; i < m; i++)
      {
        Site scaled = solvent[i];          /* dont modify the input */
        scaled.charge *= scale_coul_short; /* modify a copy */

        field (BHD, scaled, n, solute, ljc, us[i]);
      }

  /*
    1.   Put  the  solute  charge density  (positive  cores,  negative
    electrons) into Vec uc_rho.

    The density of the solute as a superposition of gaussian cores and
    distributed (electron) charge density. In  the MM case there is no
    electron density, and the core charges correspond to the effective
    (small) charges of atoms in  the solute molecule.  In the QM case,
    however,  the  core charges  (partially)  compensate the  electron
    density  and may  become large  (e.g.  +6  or +4  for  oxygen atom
    depending on the ECP used).

    The  overall factor for  the Coulomb  potential derived  from this
    density will be 1.0, that is idependent of the solvent charge.
  */

  /*
    a)  Core density  is a  superposition  of Gaussians.   For such  a
    density   the    corresponding   Coulomb   potential    is   known
    analytically. So we  have a choice to either  treat Gaussian cores
    as a most general smoothly  varying density and plug that into the
    Poisson equation or use the analytic properties. Flip this flag to
    choose:
  */
  const bool analytic_form = true;
  {
    real q[n];              /* solute charges */
    real x[n][3];           /* solute coordinates */

    /* Copy them from struct Site: */
    for (int i = 0; i < n; i++)
      {
        q[i] = solute[i].charge;
        FOR_DIM
          x[i][dim] = solute[i].x[dim];
      }

    /*
      Compute  the charge  density  of gaussian  broadened cores.   We
      primarily need the potential, but charge density is also used to
      compute some observables:
    */
    if (uc_rho)
      cores (BHD, n, q, x, G, uc_rho);

    /*
      Real-space representation of  the long-range potential is needed
      to compute the chemical potential.  One should not assume that a
      Fourier transform  of uc_fft will  deliver something meaningfull
      for charged species:
    */
    if (uc)
      coulomb_long (BHD, n, q, x, G, uc);

    if (uc_fft && analytic_form)
      coulomb_long_fft (BHD, n, q, x, G, uc_fft);
    /* else: obtain it from uc_rho, see below ... */
  }

  /* In this case we alredy got both, density and FFT of the
     potential: */
  if (analytic_form)
    {
      /* FIXME: need to combine with the field of electrons */
      assert (density == NULL);
      return;
    }

  /* b) Electron density.  Only if not NULL, add the charge density of
     electrons: */
  if (density)
    {
      assert (uc_rho != NULL);
      local Vec rho_elec = vec_duplicate (uc_rho);

      /*
        Use  density() to  compute  total electron  charge density  at
        every point  of the local grid  portion and put  that into Vec
        rho_elec:
      */
      grid_map (BHD->da, BHD->PD, density, rho_elec);

      /* Electrons are negative: */
      VecAXPY (uc_rho, -1.0, rho_elec);

      vec_destroy (&rho_elec);
    }

  /*
    This branch  was rarely tested.   Replace the computed  density by
    the one from file. We  use bgy3d_vec_read() here so the file needs
    to   be   in   the    PETSC   binary   Vec   format.    See   also
    bgy3d_vec_read_gpaw().
  */
  {
    char filename[260];         /* electron density file */

    if (bgy3d_getopt_string ("load-charge", sizeof filename, filename))
      {
        assert (uc_rho != NULL);
        bgy3d_vec_read (filename, uc_rho);
      }
  }


  if (uc_rho)                   /* debug prints only */
    {
      PetscScalar sum;
      real dV = volume_element (BHD->PD);
      VecSum (uc_rho, &sum);
      PRINTF ("integrated charge = %f (should be close to zero)\n",
                   sum * dV);
    }

  /*
    2.  Solve the Poisson equation, Δu = -4πρ/ε₀.

    Note  that the  nuclei  were gaussian-smeared  for the  long-range
    electrostatics.   Therefore  the  Vec  uc_rho  is,  formally,  not
    representing the  whole of the  solute charge density,  rather its
    "diffuse" or "smeared" part.

    Note that we only get FFT coulomb potential from Poisson solver
   */
  if (uc_fft)
    {
      assert (uc_rho != NULL);
      bgy3d_poisson (BHD, uc_fft, uc_rho, -4 * M_PI * EPSILON0INV);
    }
}
