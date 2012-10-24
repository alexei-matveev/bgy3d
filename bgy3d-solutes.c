/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O_solutes.c,v 1.3 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

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
#include "bgy3d-solvents.h"
#include "bgy3d-pure.h"         /* Coulomb_short() */
#include "bgy3d-getopt.h"
#include "bgy3d-solutes.h"

/* Solute is  isomorphic to an  array of sites.  Consider  handling it
   like that in the code.   Structs with flexible array members may be
   confusing. Such struct is convenient for literal data, though: */
typedef struct Solute {
  const char *name;             /* human readable name */
  int n;                        /* number of sites */
  Site sites[];                 /* site descriptions */
} Solute;

static void poisson (const State *BHD, Vec uc, Vec rho, real q);

/*
 * This function  obeys the callback interface assumed  in field(). It
 * is supposed to  get (1) parameters of the solvent  site such as its
 * location and force  field parameters, and (2) a  description of the
 * solute and  return a real number  such as an  interaction energy or
 * the charge density:
 */
static real ljc (const Site *A, int n, const Site S[n]);

/*
 * This function expects a callback obeying the above interface as one
 * of the arguments:
 */
static void field (DA da, const ProblemData *PD,
                   Site A, int n, const Site S[n],
                   real fact,
                   real (*f)(const Site *A, int n, const Site S[n]),
                   Vec v);

/* This is  another type of  (elemental) callbacks that take  array of
   coordinates and return an array of respective values: */
static void gf_density (int m, const real x[m][3], real rho[m],
                        int n, const Site S[n]); /* extra */

/* Callback here  is that of  gf_density() closure over all  but three
   leading arguments: */
static void grid_map (DA da, const ProblemData *PD,
                      void (*f)(int m, const real x[m][3], real fx[m]),
                      Vec v);

static void read_charge_density (DA da, const ProblemData *PD,
				const char *filename, real fact, Vec v);

// FIXME: maybe #include "solutes.h" instead?

/*********************************/
/* HCl */
/*********************************/

static const Solute HydrogenChloride =
  {"hydrogen chloride", 2,
   {{"H", {0.6285, 0.0, 0.0}, 2.735, 0.03971, 0.2},
    {"Cl", {-0.6285, 0.0, 0.0}, 3.353, 0.51434, -0.2}}};

/*********************************/
/* CS2 */
/*********************************/

static const Solute CarbonDisulfide =
  {"carbon disulfide", 3,
   {{"C", {0.0, 0.0, 0.0}, 3.2, 0.10128, -0.308},
    {"S1", {-1.56, 0.0, 0.0}, 3.52, 0.395, 0.154},
    {"S2", {1.56, 0.0, 0.0}, 3.52, 0.395, 0.154}}};

/*********************************/
/* Water */
/*********************************/

static const Solute Water =
  {"water", 3,
   {{"O", {-0.2929, 0.0, 0.0}, 3.1506, 0.1521, -0.834},
    {"OH", {0.2929, 0.757, 0.0}, 0.4, 0.046, 0.417},
    {"OH", {0.2929, -0.757, 0.0}, 0.4, 0.046, 0.417}}};

/*********************************/
/* Methanol */
/*********************************/

static const Solute Methanol =
  {"methanol", 6,
   {{"C", {-0.748, -0.015, 0.024}, 3.5, 0.066, 0.145},
    {"HC1", {-1.293, -0.202, -0.901}, 2.5, 0.03, 0.04},
    {"HC2", {-1.263, 0.754, 0.6}, 2.5, 0.03, 0.04},
    {"HC3", {-0.699, -0.934, 0.609}, 2.5, 0.03, 0.04},
    {"O", {0.558, 0.42, -0.278}, 3.12, 0.17, -0.683},
    {"OH", {0.716, 1.404, 0.137}, 0.4, 0.04, 0.418}}};

/* BUTANOIC ACID */
/* H1 sigma and epsilon adopted */

static const Solute ButanoicAcid =
  {"butanoic acid", 14,
   {{"C1", {1.422, -0.017, 0.0}, 3.75, 0.105, 0.52},
    {"O1", {1.422, 1.353, 0.0}, 2.96, 0.21, -0.44},
    {"O2", {2.643, -0.722, 0.0}, 3.0, 0.17, -0.53},
    {"C2", {0.1, -0.78, 0.0}, 3.5, 0.066, -0.12},
    {"C3", {-1.06, 0.212, 0.0}, 3.5, 0.066, -0.12},
    {"C4", {-2.381, -0.551, 0.0}, 3.5, 0.066, -0.18},
    {"OH", {3.21, -0.461, 0.882}, 3.4, 0.046, 0.45},
    {"H2", {0.043, -1.407, 0.89}, 2.5, 0.03, 0.06},
    {"H3", {0.043, -1.407, -0.89}, 2.5, 0.03, 0.06},
    {"H4", {-1.002, 0.838, -0.89}, 2.5, 0.03, 0.06},
    {"H5", {-1.002, 0.838, 0.89}, 2.5, 0.03, 0.06},
    {"H6", {-2.439, -1.178, 0.89}, 2.5, 0.03, 0.06},
    {"H7", {-2.439, -1.178, -0.89}, 2.5, 0.03, 0.06},
    {"H8", {-3.21, 0.157, 0.0}, 2.5, 0.03, 0.06}}};

/*********************************/
/* Hexane */
/*********************************/

static const Solute Hexane =
  {"hexane", 20,
   {{"C", {1.709, -2.812, 0.0}, 3.5, 0.066, -0.18},
    {"C", {1.684, -1.278, 0.0}, 3.5, 0.066, -0.12},
    {"C", {0.245, -0.753, 0.0}, 3.5, 0.066, -0.12},
    {"C", {0.241, 0.779, 0.0}, 3.5, 0.066, -0.12},
    {"C", {-1.198, 1.304, 0.0}, 3.5, 0.066, -0.12},
    {"C", {-1.206, 2.834, 0.0}, 3.5, 0.066, -0.18},
    {"H", {2.236, -3.164, 0.887}, 2.5, 0.03, 0.06},
    {"H", {2.232, -3.164, -0.89}, 2.5, 0.03, 0.06},
    {"H", {0.691, -3.204, 0.003}, 2.5, 0.03, 0.06},
    {"H", {2.202, -0.914, -0.888}, 2.5, 0.03, 0.06},
    {"H", {2.201, -0.914, 0.89}, 2.5, 0.03, 0.06},
    {"H", {-0.273, -1.115, 0.889}, 2.5, 0.03, 0.06},
    {"H", {-0.272, -1.115, -0.89}, 2.5, 0.03, 0.06},
    {"H", {0.757, 1.142, -0.89}, 2.5, 0.03, 0.06},
    {"H", {0.757, 1.141, 0.89}, 2.5, 0.03, 0.06},
    {"H", {-1.716, 0.944, 0.89}, 2.5, 0.03, 0.06},
    {"H", {-1.716, 0.944, -0.89}, 2.5, 0.03, 0.06},
    {"H", {-0.696, 3.204, -0.89}, 2.5, 0.03, 0.06},
    {"H", {-0.696, 3.204, 0.89}, 2.5, 0.03, 0.06},
    {"H", {-2.236, 3.19, 0.0}, 2.5, 0.03, 0.06}}};

static const Solute *solutes[] = {&HydrogenChloride, /* 0 */
                                  &CarbonDisulfide,  /* 1 */
                                  &Water,            /* 2 */
                                  &Methanol,         /* 3 */
                                  &ButanoicAcid,     /* 4 */
                                  &Hexane};          /* 5 */

/* Get solute  sites by solute name.  Functions that do  the real work
   operate on array of sites: */
void bgy3d_solute_get (const char *name, int *n, const Site **sites)
{
  /* Number of known solutes: */
  const size_t len = sizeof (solutes) / sizeof (*solutes);

  /* Check whether we have the solute name in solutes[]: */
  size_t i = 0;
  while (i < len && strcmp (name, solutes[i]->name))
    i++;

  /* Could not find the solute name: */
  assert (i < len);

  /* Return solute sites and their count: */
  *n = solutes[i]->n;
  *sites = solutes[i]->sites;
}

/*
 * Create initial solute data.
 *
 * See (5.106) and (5.08) in the  thesis.  Fill us[0] and us[1], for H
 * and O in that order  with short-range, VM_LJ + VM_coulomb_short and
 * uc  with long range  potential VM_coulomb_long  (not scaled  by the
 * charges of the solvent sites).
 *
 * The function  pointer density() passed  to bgy3d_solute_field(), if
 * not NULL,  may be used  to compute the distributed  charge density,
 * e.g. originating from electrons.  FIXME: be more specific about the
 * sign, is it a charge density or the electron density?
 */

void bgy3d_solute_field (const State *BHD,
                         int m, const Site solvent[m], /* m == 2 */
                         Vec us[m], Vec uc, /* intent(out) */
                         int n, const Site solute[n], /* n arbitrary */
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         real damp, real damp_LJ)
{
  PetscPrintf (PETSC_COMM_WORLD,
               "Computing solute data with damping factors %f, %f\n",
               damp, damp_LJ);


  /*
    Calculate FF potential for all solvent sites.

    FIXME:  scaling the  epsilon of  the  solvent site  by factor  X^2
    scales the interaction with the  solute by factor X. Why not using
    this fact instead of handling damp/damp_LJ separately?  Similarly,
    scaling the  solvent site  charge by a  factor scales  the Coulomb
    interaction.  At  least in the  two test examples the  two factors
    are  identical.   Initial  version  of  the code  scaled  LJ-  and
    short-range Coulomb interaction  using two different factors.  The
    new code uses just one, that is why the assertion:
  */
  real factor = damp;
  assert (damp == damp_LJ);
  assert (factor >= 0.0);

  /*
   * Fill the force field interaction of H and O (or other) sites with
   * the solute into the respective arrays.
   *
   * We  supply ljc()  as  a  callback function  that  is supposed  to
   * compute the  interaction of  a charged LJ  solvent site  with the
   * solute.
   */
#ifndef QM
  for (int i = 0; i < m; i++)
    field (BHD->da, BHD->PD, solvent[i], n, solute, factor, ljc,
           us[i]);
#else
  /* At  this  place the  (short  range)  Coulomb  interaction of  the
    solvent  site   with  the  solute  was   deliberately  omitted  by
    specifying zero charge of the solvent site. This effectively makes
    a point  charge (Coulomb  short + Coulomb  long) to  a distributed
    Gaussian (Coulomb long only). */

  for (int i = 0; i < m; i++)
    {
      Site neutral = solvent[i]; /* dont modify the global variable */
      neutral.charge = 0.0;      /* modify a copy */

      field (BHD->da, BHD->PD,
             neutral, n, solute, factor, ljc,
             us[i]);
    }
#endif

  /*
   * Compute the charge density  of the solute.  The callback function
   * gf_density() sums  (net or nuclear) charge  distribution for each
   * solute site. The overall factor for the Coulomb potential derived
   * from this density is 1.0 (idependent of the solvent charge).
   *
   * Vec uc  will first hold the  solute charge density  and later its
   * Coulomb field.
   *
   * 1.  Put  the  solute  charge  density  (positive  core,  negative
   * electrons) into Vec uc:
   */

  char filename[260];           /* electron density file */

  if (bgy3d_getopt_string("--load-charge", filename, sizeof (filename)))
    {
      read_charge_density (BHD->da, BHD->PD, filename, 1.0, uc);
    }
  else
    {
      /* This function computes the density  of the solute as a sum of
         gaussian functions and  distributed (electron) charge density
         at an  array of points.  FIXME: nested  closure function, GCC
         extension: */
      void f (int m, const real x[m][3], real rho[m])
      {
        /* Bind solute description  from the enclosing scope.  Compute
           the  (positive)   charge  density  of  (gaussian-broadened)
           cores: */
        gf_density (m, x, rho, n, solute);

        /* If not NULL, add charge density of electrons: */
        if (density)
          {
            /* Electron density goes here: */
            real rho1[m];       /* MEMORY: huge array here! */

            /* This   computes   the   (unsigned)   electron   density
               distributed in space by QM rules: */
            density (m, x, rho1);

            for (int i = 0; i < m; i++)
              rho[i] -= rho1[i]; /* electrons are negative */
          }
      }

      /* Use f()  to compute total core  & electron at  every point of
         the local grid portion and put that into Vec uc: */
      grid_map (BHD->da, BHD->PD, f, uc);
    }

  if (1)                        /* debug prints only */
    {
      PetscScalar sum;
      real dV = BHD->PD->h[0] * BHD->PD->h[1] * BHD->PD->h[2];
      VecSum (uc, &sum);
      PetscPrintf (PETSC_COMM_WORLD,
                   "integrated charge = %f (should be close to zero)\n",
                   sum * dV);
    }

  /*
   * 2. Solve  the Poisson equation "in-place" by  specifying the same
   * Vec uc as input and output:
   */
  poisson (BHD, uc, uc, 1.0 * damp); /* WARNING: argument aliasing here! */
}

/*
 * Calculate a  real field "f"  for the solvent site  characterized by
 * (sigma,  epsilon, charge) in the presence  of the solute  S with an
 * overall factor "fact" at every point (x, y, z) of the local grid.
 *
 * The function f(x, y, z, eps, sig,  chg, S) can be ljc() or rho() as
 * two examples.
 *
 * Site  A  is intent(in)  but  has to  be  passed  by value  (copying
 * involved) at the  moment.  We are modifying coordinates  of a local
 * copy  and pass a  reference to  that copy  further to  the callback
 * function *f.  The coordinates of a site as passed to this functions
 * are ignored.
 *
 * Vector "v" is the intent(out) argument.
 */
static void field (DA da, const ProblemData *PD,
                   Site A, int n, const Site S[n],
                   real fact,
                   real (*f)(const Site *A, int n, const Site S[n]),
                   Vec v)
{
  PetscScalar ***vec;
  real h[3];
  int i0, j0, k0;
  int ni, nj, nk;

  /*
   * FIXME: do we  really assume that intervals for x-,  y- and z- are
   * the same? This basically means the  corner of the unit cell is at
   * (offset, offset, offset):
   */
  real offset = PD->interval[0];

  FOR_DIM
    h[dim] = PD->h[dim];

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  DAVecGetArray (da, v, &vec);

  /* loop over local portion of grid */
  for (int k = k0; k < k0 + nk; k++)
    {
      A.x[2] = k * h[2] + offset;

      for (int j = j0; j < j0 + nj; j++)
        {
          A.x[1] = j * h[1] + offset;

          for (int i = i0; i < i0 + ni; i++)
            {
              A.x[0] = i * h[0] + offset;

              /*
               * Compute the field  f at (x, y, z) <->  (i, j, k) e.g.
               * by summing  (LJ) contributions from  all solute sites
               * at that grid point:
               */
              vec[k][j][i] = fact * f (&A, n, S);
            }
        }
    }
  DAVecRestoreArray (da, v, &vec);
}

static void grid_map (DA da, const ProblemData *PD,
                      void (*f)(int m, const real x[m][3], real fx[m]),
                      Vec v)
{
  int i0, j0, k0;
  int ni, nj, nk;

  /* Get local portion of the grid */
  DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

  /* MEMORY: huge arrays here: */
  int m = ni * nj * nk;
  real x[m][3], fx[m];

  /* Get coordinates of the local grid portion: */
  {
    int ijk = 0;
    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          {
            /* Coordinates  (x, y,  z) <->  (i, j,  k).  FIXME:  do we
               really assume that intervals for  x-, y- and z- are the
               same? This basically means  the corner of the unit cell
               is at (offset, offset, offset): */
            x[ijk][0] = i * PD->h[0] + PD->interval[0];
            x[ijk][1] = j * PD->h[1] + PD->interval[0];
            x[ijk][2] = k * PD->h[2] + PD->interval[0];
            ijk++;
          }
    assert (ijk == m);
  }

  /* Let the function f() compute  the field/density at our portion of
     the grid. FIXME: cast here is to silence the const-warning: */
  f (m, (const real (*)[3]) x, fx);

  /* Copy contents of fx[] to the output vector: */
  {
    PetscScalar ***vec;
    DAVecGetArray (da, v, &vec);
    int ijk = 0;
    for (int k = k0; k < k0 + nk; k++)
      for (int j = j0; j < j0 + nj; j++)
        for (int i = i0; i < i0 + ni; i++)
          {
            vec[k][j][i] = fx[ijk];
            ijk++;
          }
    assert (ijk == m);
    DAVecRestoreArray (da, v, &vec);
  }
}

/*
 * Interaction of a charged LJ site (sigma, epsilon, charge) at (x, y,
 * z) with the solute S:
 */
static real ljc (const Site *A, int n, const Site S[n])
{
  /* Sum force field contribution from all solute sites: */
  real field = 0.0;

  for (int site = 0; site < n; site++)
    {

      /* Interaction parameters for a pair of LJ sites: */
      real e2 = sqrt (A->epsilon * S[site].epsilon);
      real s2 = 0.5 * (A->sigma + S[site].sigma);

      /* Distance from a grid point to this site: */
      real r_s = sqrt (SQR(A->x[0] - S[site].x[0]) +
                       SQR(A->x[1] - S[site].x[1]) +
                       SQR(A->x[2] - S[site].x[2]));

      /* 1. Lennard-Jones */
      field += Lennard_Jones (r_s, e2, s2);

      /* 2.  Coulomb, short  range part.   For historical  reasons the
         overall scaling  factor, the  product of solvent-  and solute
         site charges, is handled by the function itself: */
      field += Coulomb_short (r_s, A->charge * S[site].charge);
    }

  return field;
}

/*
 * Gaussian functions (gf) represent  each solute site charge, compute
 * the  total  charge density  at  every  point  x[]. Compare  to  the
 * function rho() above that does the same for one point at a time.
 *
 * Each gaussian is evaluated as:
 *
 *   ρ(r) = q * [G / √π]³ * exp[-G² * (r - x₀)²]
 */
static void gf_density (int m, const real x[m][3], /* coordinates */
                        real rho[m],            /* output densities */
                        int n, const Site S[n]) /* solute description */
{
  /* G is predefind in bgy3d-solvents.h FIXME: make the gaussian width
     a property of the (solute) site  in the same way as the charge of
     the site. */
  real prefac = pow(G / sqrt(M_PI), 3.0);

  for (int i = 0; i < m; i++)
    {
      /* Sum Gaussian contributions from all solute sites: */
      real ro = 0.0;
      for (int j = 0; j < n; j++)
        {

          /* Square of the distance from a grid point to this site: */
          real r2 = SQR (x[i][0] - S[j].x[0]) +
                    SQR (x[i][1] - S[j].x[1]) +
                    SQR (x[i][2] - S[j].x[2]);

          /* Gaussian distribution,  note that G  is not a  width, but
             rather an inverse of it: */
          ro += prefac * S[j].charge * exp(- G * G * r2);
        }
      rho[i] = ro;
    }
}


/*
  Solve  Poisson  Equation  in  Fourier space  and  get  elestrostatic
  potential by inverse FFT.

  Vec uc is intent(out).
  Vec rho is intent(in).
  real q is the overall factor.

  As a  matter of  fact, it  appears that one  could provide  the same
  factual parameter  for rho and  uc to effectively solve  the Poisson
  equation "in place".

  Except of  temporary allocation of a  complex Vec does  not have any
  side effect.
*/
void poisson (const State *BHD, Vec uc, Vec rho, real q)
{
  const real *interval = BHD->PD->interval; /* [2] */
  const real *h = BHD->PD->h;   /* h[3] */
  const int *N = BHD->PD->N;    /* N[3] */

  const real L = interval[1] - interval[0];
  const real h3 = h[0] * h[1] * h[2];

  /* Scratch complex vector: */
  Vec work;
  DACreateGlobalVector (BHD->dc, &work);

  /* Get FFT of  rho: rho(i, j, k) -> fft_rho(kx,  ky, kz) placed into
     complex work: */
  MatMult (BHD->fft_mat, rho, work);

  /*
    Solving Poisson Equation (SI units) with FFT and IFFT:

        - ΔU (x, y, z) = (1 / ε₀) ρ(x, y, z)
            c

    because of x = i h, y = j h, and z = k h, with grid spacing h =
    L/n:

        - n² / L²  Δuc(i, j, k) = (1 / ε₀) ρ(i, j, k)

    FFT (see FFTW manual "What FFTW Really Computes"):

    fft_uc(kx, ky, kz) = 1 / [4 π²  ε₀ k²  / L² ] fft_rho(kx, ky, kz)

    with

        k² = kx² + ky² + kz²

    IFFT (see FFTW manual "What FFTW Really Computes"):

    because: IFFT(fft_uc(kx, ky, kz)) = n³ * uc(i, j, k)

    uc(i, j, k) = h³ / L³  * IFFT(fft_uc(kx, ky, kz))
  */

  /* EPSILON0INV = 1 / 4 π ε₀: */
  const real scale = q * EPSILON0INV / M_PI * h3 / (L * L * L);

  /* Loop over local portion of the k-grid */
  {
    int x[3], n[3], i[3], ic[3];
    DAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

    struct {PetscScalar re, im;} ***work_;
    DAVecGetArray (BHD->dc, work, &work_);

    for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
      for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
        for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
          {
            /* FIXME: what is we  change the complex vectors to remove
               the redundnacy? */
            FOR_DIM
              {
                if (i[dim] <= N[dim] / 2)
                  ic[dim] = i[dim];
                else
                  ic[dim] = i[dim] - N[dim];
              }

            if (ic[0] == 0 && ic[1] == 0 && ic[2] == 0)
              {
                /* The gamma point, k = 0, we cannot divide by 0: */
                work_[i[2]][i[1]][i[0]].re = 0.0;
                work_[i[2]][i[1]][i[0]].im = 0.0;
              }
            else
              {
                const real k2 = (SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0])) / SQR(L);

                const real fac = scale / k2;

                /* Here we compute in place: uc(kx, ky, kz) := scale *
                   rho(kx, ky, kz) / k^2 */
                work_[i[2]][i[1]][i[0]].re *= fac;
                work_[i[2]][i[1]][i[0]].im *= fac;
              }
          }
    DAVecRestoreArray (BHD->dc, work, &work_);
  }

  /* uc := IFFT(uc(kx, ky, kz)) */
  MatMultTranspose (BHD->fft_mat, work, uc);

  VecDestroy (work);
}

/* Read density file generated by  GPAW calculation need BOHR to scale
   the length values back. */
static void read_charge_density (DA da, const ProblemData *PD,
                                 const char *filename, real fact, Vec v)
{
  PetscScalar ***vec;
  char line_buffer[BUFSIZ];
  int AtomNum;                  /* atom numbers in the molecule */
  real corner[3], dx[3], dy[3], dz[3];
  int GridNum[3]; /* Grid numers in each direction: [0]:x, [1]:y, [2]:z */
  real h[3];
  FILE *fp;
  int i0, j0, k0;
  int ni, nj, nk;


  fp = fopen(filename, "r");
  if (fp == NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Can not open file %s. \n", filename);
      exit(1);
    }

  PetscPrintf(PETSC_COMM_WORLD, "Reading data from %s. \n", filename);
  /* Skip the first two lines */
  for (int i = 0; i < 2; i++)
    {
      char *p = fgets (line_buffer, sizeof(line_buffer), fp);
      assert (p != NULL);
    }

  /* Atom numers and corner shifts */
  {
    int n = fscanf (fp, "%d %lf %lf %lf", &AtomNum, &corner[0], &corner[1], &corner[2]);
    assert (n == 4);
  }

  /* Grid numbers and spaces */
  FOR_DIM
    {
      int n = fscanf (fp, "%d %lf %lf %lf", &GridNum[dim], &dx[dim], &dy[dim], &dz[dim]);
      assert (n == 4);
    }

  /* Scale the grid space */
  dx[0] *= BOHR;
  dy[1] *= BOHR;
  dz[2] *= BOHR;

  /* Allocate memory */
  int electron[AtomNum];  /* electrons of each atom */
  real zero[AtomNum];     /* = 0.0 as in ase.io.cube.write_cube, don't
                             know the meaning */
  real x[AtomNum], y[AtomNum], z[AtomNum];

  for (int i = 0; i < AtomNum; i++)
    {
      int n = fscanf (fp, "%d %lf %lf %lf %lf", &electron[i], &zero[i], &x[i], &y[i], &z[i]);
      assert (n == 5);
      x[i] *= BOHR;
      y[i] *= BOHR;
      z[i] *= BOHR;
    }

  DAGetCorners(da, &i0, &j0, &k0, &ni, &nj, &nk);

  FOR_DIM
    h[dim] = PD->h[dim];


  /* FIXME: will need interpolation if grid not match */
  if ( GridNum[0] * GridNum[1] * GridNum[2] != ni * nj * nk)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Grid size not match!\n");
      exit(1);
    }
  else if ( fabs(dx[0] - h[0]) >= 0.001 || fabs(dy[1] - h[1]) >= 0.001 || fabs(dz[2] - h[2]) >= 0.001)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Grid space not match!\n");
      exit(1);
    }

  DAVecGetArray(da, v, &vec);

  /* electron density scaled by BOHR^3 in python script */
  real invB3 = 1. / BOHR / BOHR / BOHR;
  for (int i = i0; i < i0 + ni; i++)
    for (int j = j0; j < j0 + nj; j++)
      for (int k = k0; k < k0 + nk; k++)
        {
          int n = fscanf (fp, "%lf", &vec[i][j][k]);
          assert (n == 1);
          vec[i][j][k] *= fact * invB3;
        }

  fclose(fp);
  DAVecRestoreArray(da, v, &vec);
}

