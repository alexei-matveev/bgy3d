/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/*
 * This file  was initially named bgy3dH2O_solutes.c, as  you may have
 * guessed from  the Id  string above. There  is not much  specific to
 * water in  this file.  Rather it is  general for all  2-site solvent
 * model. So it was later renamed.
 *
 * The  funciton  bgy3d_solute_field()  initiates the  computation  of
 * (short  range)  site-specifc force  field  and  long range  Coulomb
 * potential (common  for all sites) for arbitrary  number of solvent-
 * and solute sites.
 */

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-poisson.h"      /* bgy3d_poisson() */
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* G_COULOMB_INVERSE_RANGE */
#include "bgy3d-force.h"        /* lennard_jones_coulomb_short() */
#include "bgy3d-vec.h"          /* bgy3d_vec_read() */


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
      real e2 = sqrt (A->epsilon * B->epsilon);
      real s2 = 0.5 * (A->sigma + B->sigma);
      real q2 = A->charge * B->charge;

      /* Distance from a grid point to this site: */
      real r_s = sqrt (SQR(A->x[0] - B->x[0]) +
                       SQR(A->x[1] - B->x[1]) +
                       SQR(A->x[2] - B->x[2]));

      /* Lennard-Jones + Coulomb, short range part: */
      field += lennard_jones_coulomb_short (r_s, s2, e2, G, q2);
    }

  return field;
}


/*
  Gaussian functions  represent each  solute site charge,  compute the
  total charge density at every point x[i].

  Each gaussian is evaluated as:

    ρ(r) = q * [G / √π]³ * exp[-G² * (r - x₀)²]
*/
static void
cores (const State *BHD,
       int n, const real q[n], /* const */ real r[n][3], real G,
       Vec rho)                /* out, real, center */
{
  /*
    FIXME: make the gaussian width  a property of the (solute) site in
    the same way as the charge of the site.
  */
  const real prefac = pow (G / sqrt (M_PI), 3.0);

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
        sum += q[i] * exp(- G * G * r2);
      }
    return prefac * sum;
  }
  vec_rmap3 (BHD, f3, rho);
}


/*
  Tabulate ρ(k) for ρ(r) being  a superposition of δ-like functions at
  r[i]  with weights q[i].  Centered at  the corner.  Normalization is
  such that

    ρ(k=0) = Σ q[i]
              i

  FIXME: if this ever  becomes performance critical consider using the
  properties of the exponential and maybe constant grid spacing.

     a + b    a    b
    e      = e  * e
*/
static void
deltas_fft (const State *BHD,
            int m, const real q[m], real r[m][3], /* in */
            Vec rho_fft)        /* out, complex, corner */
{
  complex f3 (const real k[3])
  {
    complex sum = 0.0;
    for (int p = 0; p < m; p++)
      sum += q[p] * cexp (-I * (k[0] * r[p][0] + k[1] * r[p][1] + k[2] * r[p][2]));
    return sum;
  }
  vec_kmap3 (BHD, f3, rho_fft);
}


/*
  Return ρ(k) for ρ(r) being  a superposition of gaussian at r[i] with
  weights q[i].  Centered at the corner.  Normalization is such that

    ρ(k=0) = Σ q[i]
              i

  Each gaussian is evaluated as:

    ρ(k) = q * exp[-k² / 4G²]

  The phase factor accounts for the translation.
*/
static void
cores_fft (const State *BHD,
           int n, const real q[n], /* const */ real r[n][3], real G,
           Vec rho_fft)            /* out, complex, corner */
{
  /* Compute the k-representation of the position density: */
  deltas_fft (BHD, n, q, r, rho_fft);

  /* K-representation of the (gaussian) core shape f(r) goes here: */
  local Vec f_fft = vec_duplicate (rho_fft);

  /*
    This is  the k-representation of  a gaussian unit charge  of width
    1/G:
                    2    2
      f(k) = exp (-k / 4G )

    We will apply it as a convolution kernel to the point-density.
   */
  complex f (real k)
  {
    return exp (- SQR (k) / (4 * SQR (G)));
  }
  vec_kmap (BHD, f, f_fft);

  /* Pointwise multiplication of complex numbers: */
  complex pure mul (complex x, complex y)
  {
    return x * y;
  }
  vec_fft_map2 (rho_fft, mul, f_fft, rho_fft); /* aliasing! */

  /*
    FIXME: note how we constructed a temporary Vec f_fft just to scale
    rho_fft  by  a  known  function  of  k ---  the  Vec  rho_fft  has
    realtively simple closed form.
  */
  vec_destroy (&f_fft);
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
                    Vec us[m],                    /* out */
                    Vec uc_fft,                   /* out, complex */
                    Vec uc_rho,                   /* out, optional */
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
  for (int i = 0; i < m; i++)
    {
      Site scaled = solvent[i];          /* dont modify the input */
      scaled.charge *= scale_coul_short; /* modify a copy */

      field (BHD, scaled, n, solute, ljc, us[i]);
    }

  /* Early return if no Coulomb is requested: */
  if (uc_fft == NULL && uc_rho == NULL)
    return;

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

  /* a) Core density: */
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

    /* Compute the charge density of gaussian broadened cores: */
    cores (BHD, n, q, x, G, uc_rho);
  }

  /* b) Electron density.  Only if not NULL, add the charge density of
     electrons: */
  if (density)
    {
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
      bgy3d_vec_read (filename, uc_rho);
  }


  if (true)                     /* debug prints only */
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
  bgy3d_poisson (BHD, uc_fft, uc_rho, -4 * M_PI * EPSILON0INV);
}
