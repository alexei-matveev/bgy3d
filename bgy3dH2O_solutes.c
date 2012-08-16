/*==========================================================*/
/*  $Id: bgy3dH2O_solutes.c,v 1.3 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d_SolventParameters.h"

typedef struct Site {
    char name[5];            /* atom types. What are they used for? */
    real x[3];               /* coordinates */
    real sigma;              /* sigma for LJ */
    real epsilon;            /* epsilon for LJ */
    real charge;             /* charge */
} Site;

/* Solute is  isomorphic to an  array of sites.  Consider  handling it
   like that in the code.   Structs with flexible array members may be
   confusing. Such struct is convenient for literal data, though: */
typedef struct Solute {
    int n;                      /* number of sites */
    Site sites[];               /* site descriptions */
} Solute;

static void poisson (BGY3dH2OData BHD, Vec uc, Vec rho, real q);
static void solute_field (BGY3dH2OData BHD, const Site S[], int nsites, real damp, real damp_LJ);

/*
 * These two functions  obey the same interface. They  are supposed to
 * get (1)  parameters of  the solvent site  such as its  location and
 * force field  parameters, and  (2) a description  of the  solute and
 * return a  real number such as  an interaction energy  or the charge
 * density:
 */
static real ljc (const Site *A, const Site S[], int nsites);
static real rho (const Site *A, const Site S[], int nsites);

/*
 * This function expects a callback obeying the above interface as one
 * of the arguments:
 */
static void field (DA da, const ProblemData *PD,
                   Site A, const Site S[], int nsites,
                   real fact,
                   real (*f)(const Site *A, const Site S[], int nsites),
                   Vec v);

// FIXME: maybe #include "solutes.h" instead?

/*********************************/
/* Water */
/*********************************/

static const Solute Water =
  {3, {{"O", {-0.2929, 0.0, 0.0}, 3.1506, 0.1521, -0.834},
       {"OH", {0.2929, 0.757, 0.0}, 0.4, 0.046, 0.417},
       {"OH", {0.2929, -0.757, 0.0}, 0.4, 0.046, 0.417}}};

/*********************************/
/* CS2 */
/*********************************/

static const Solute CarbonDisulfide =
  {3, {{"C", {0.0, 0.0, 0.0}, 3.2, 0.10128, -0.308},
       {"S1", {-1.56, 0.0, 0.0}, 3.52, 0.395, 0.154},
       {"S2", {1.56, 0.0, 0.0}, 3.52, 0.395, 0.154}}};

/*********************************/
/* HCl */
/*********************************/

static const Solute HydrogenChloride =
  {2, {{"H", {0.6285, 0.0, 0.0}, 2.735, 0.03971, 0.2},
       {"Cl", {-0.6285, 0.0, 0.0}, 3.353, 0.51434, -0.2}}};

/*********************************/
/* Methanol */
/*********************************/

static const Solute Methanol =
  {6, {{"C", {-0.748, -0.015, 0.024}, 3.5, 0.066, 0.145},
       {"HC1", {-1.293, -0.202, -0.901}, 2.5, 0.03, 0.04},
       {"HC2", {-1.263, 0.754, 0.6}, 2.5, 0.03, 0.04},
       {"HC3", {-0.699, -0.934, 0.609}, 2.5, 0.03, 0.04},
       {"O", {0.558, 0.42, -0.278}, 3.12, 0.17, -0.683},
       {"OH", {0.716, 1.404, 0.137}, 0.4, 0.04, 0.418}}};

/*********************************/
/* Hexane */
/*********************************/

static const Solute Hexane =
  {20, {{"C", {1.709, -2.812, 0.0}, 3.5, 0.066, -0.18},
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

/* BUTANOIC ACID */
/* H1 sigma and epsilon adopted */

static const Solute ButanoicAcid =
  {14, {{"C1", {1.422, -0.017, 0.0}, 3.75, 0.105, 0.52},
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

/* These are  the two  solvent sites.  Coordinates  will not  be used.
   Respective parameters are #defined  elsewhere. Also do not take the
   names  of  the sites  literally.  The  same  structure is  used  to
   represent all (2-site?) solvents, such as HCl:  */
static const Site solvent[] =
  {{"h", {0.0, 0.0, 0.0}, sH, eH, qH},
   {"o", {0.0, 0.0, 0.0}, sO, eO, qO}};

static void recompute_initial_data (BGY3dH2OData BHD, const Solute *S, real damp, real damp_LJ)
{
    /* Functions that do the real work operate on array of sites: */
    solute_field (BHD, S->sites, S->n, damp, damp_LJ);
}

void RecomputeInitialSoluteData_Water(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Water.\n");
  recompute_initial_data (BHD, &Water, damp, damp_LJ);
}

void RecomputeInitialSoluteData_CS2(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is CarbonDisulfide.\n");
  recompute_initial_data (BHD, &CarbonDisulfide, damp, damp_LJ);
}

void RecomputeInitialSoluteData_HCl(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is HCl.\n");
  recompute_initial_data (BHD, &HydrogenChloride, damp, damp_LJ);
}

void RecomputeInitialSoluteData_Methanol(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Methanol.\n");
  recompute_initial_data (BHD, &Methanol, damp, damp_LJ);
}

void RecomputeInitialSoluteData_Hexane(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Hexane.\n");
  recompute_initial_data (BHD, &Hexane, damp, damp_LJ);
}

void RecomputeInitialSoluteData_ButanoicAcid(BGY3dH2OData BHD, real damp, real damp_LJ, real zpad)
{
  PetscPrintf(PETSC_COMM_WORLD,"Solute is Butanoic Acid.\n");
  recompute_initial_data (BHD, &ButanoicAcid, damp, damp_LJ);
}

/*
 * Create initial solute data.
 *
 * XXX: See (5.106) and (5.08)  in the thesis.  Return BHD->gH_ini and
 *      BHD->gO_ini (beta *  (VM_LJ + VM_coulomb_short)) and BHD->ucH,
 *      BHD->ucO (beta * VM_coulomb_long), but is beta missing here?
 */

static void solute_field (BGY3dH2OData BHD, const Site S[], int nsites, real damp, real damp_LJ)
{
  PetscPrintf(PETSC_COMM_WORLD,"Recomputing solute data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);

  /* FIXME: an alternative version (QM) did not set them, nevertheless
     the HCl case appeared to work. Are these really necessary? */
  VecSet(BHD->gHO_ini, 0.0);    /* What is it used for? */
  VecSet(BHD->ucHO, 0.0);
  FOR_DIM {
    VecSet(BHD->fH_l[dim], 0.0);
    VecSet(BHD->fO_l[dim], 0.0);
    VecSet(BHD->fHO_l[dim], 0.0); /* What is it used for? */
  }

  /*
    Calculate FF potential for all solvent sites.

    Beta  is the  (inverse) temperature.   For historical  reasons the
    solute  field  acting on  solvent  sites  is  defined having  this
    factor.

    FIXME:  scaling the  epsilon of  the  solvent site  by factor  X^2
    scales the interaction with the  solute by factor X. Why not using
    this fact instead of handling damp/damp_LJ separately?  Similarly,
    scaling the  solvent site  charge by a  factor scales  the Coulomb
    interaction.  At  least in the  two test examples the  two factors
    are  identical.   Initial  version  of  the code  scaled  LJ-  and
    short-range Coulomb interaction  using two different factors.  The
    new code uses just one, that is why the assertion:
  */
  real factor = damp * BHD->PD->beta;
  assert (damp == damp_LJ);
  assert (factor >= 0.0);

  /*
   * Fill ff-interaction  of H  and O sites  with the solute  into the
   * respective arrays.
   *
   * FIXME: Force field parameters, (eH, sH, qH) and (eO, sO, qO), for
   * H and O are #defined at some obscure place:
   *
   * We  supply ljc()  as  a  callback function  that  is supposed  to
   * compute the  interaction of  a charged LJ  solvent site  with the
   * solute.
   */
#ifndef QM
  field (BHD->da, BHD->PD, solvent[0], S, nsites, factor, ljc, BHD->gH_ini);
  field (BHD->da, BHD->PD, solvent[1], S, nsites, factor, ljc, BHD->gO_ini);
#else
  /* At  this  place the  (short  range)  Coulomb  interaction of  the
    solvent  site   with  the  solute  was   deliberately  omitted  by
    specifying zero charge of the solvent site. This effectively makes
    a point  charge (Coulomb  short + Coulomb  long) to  a distributed
    Gaussian (Coulomb long only). */

  Site neutral[2];              /* dont modify the global variable */

  for (int i = 0; i < 2; i++) {
    neutral[i] = solvent[i];
    neutral[i].charge = 0.0;    /* modify a copy */
  }

  field (BHD->da, BHD->PD, neutral[0], S, nsites, factor, ljc, BHD->gH_ini);
  field (BHD->da, BHD->PD, neutral[1], S, nsites, factor, ljc, BHD->gO_ini);
#endif

  /*
   * Compute the charge density  of the solute.  The callback function
   * rho() sums charge distribution for  each solute site and does not
   * use (epsilon,  sigma, charge) parameters of the  solvent site, so
   * that  we  provide -1.0  for  them.   The  overall factor  is  1.0
   * (idependent of the solvent charge):
   */

  Site point = {"x", {0.0, 0.0, 0.0}, -1.0, -1.0, -1.0};

  Vec v; /* Vector for solute charge density and its Coulomb field */

  /* MEMORY:  huge array  here!  FIXME:  make re-use  of pre-allocated
     vectors more transparent and get rid of this: */
  DACreateGlobalVector (BHD->da, &v);

  /* 1. Put the solute density into Vec v. Due to the inter */
  field (BHD->da, BHD->PD, point, S, nsites, 1.0, rho, v);

  /*
   * 2. Solve  the Poisson equation "in-place" by  specifying the same
   * Vec v as  input and output.  (Original version  was ouputting the
   * Coulomb potential into a pre-allocated vector BHD->v[0]).
   */
  poisson (BHD, v, v, 1.0 * damp); /* WARNING: argument aliasing here! */

  /*
   * 3. Copy  the electrostatic potential  scaled by the  solvent site
   * charges into predefined locations:
   */
  VecSet (BHD->ucH, 0.0);
  VecAXPY (BHD->ucH, qH, v);

  VecSet (BHD->ucO, 0.0);
  VecAXPY (BHD->ucO, qO, v);

  /* MEMORY: deallocate huge array here! */
  VecDestroy (v);
}

/*
 * Calculate a  real field "f"  for the solvent site  characterized by
 * (epsilon, sigma,  charge) in the presence  of the solute  S with an
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
                   Site A, const Site S[], int nsites,
                   real fact,
                   real (*f)(const Site *A, const Site S[], int nsites),
                   Vec v)
{
    PetscScalar ***vec;
    real h[3];
    int i0, j0, k0;
    int ni, nj, nk;

    /*
     * FIXME: do we really assume that intervals for x-, y- and z- are
     * the same? This  basically means the corner of  the unit cell is
     * at (offset, offset, offset):
     */
    real offset = PD->interval[0];

    FOR_DIM
        h[dim] = PD->h[dim];

    /* Get local portion of the grid */
    DAGetCorners (da, &i0, &j0, &k0, &ni, &nj, &nk);

    DAVecGetArray (da, v, &vec);

    /* loop over local portion of grid */
    for (int k = k0; k < k0 + nk; k++) {

        A.x[2] = k * h[2] + offset;

        for (int j = j0; j < j0 + nj; j++) {

            A.x[1] = j * h[1] + offset;

            for (int i = i0; i < i0 + ni; i++) {

                A.x[0] = i * h[0] + offset;

                /*
                 * Compute the field f at (x, y, z) <-> (i, j, k) e.g.
                 * by summing (LJ) contributions from all solute sites
                 * at that grid point:
                 */
                vec[k][j][i] = fact * f (&A, S, nsites);
            }
        }
    }
    DAVecRestoreArray (da, v, &vec);
}

/*
 * Interaction of a charged LJ site (epsilon, sigma, charge) at (x, y,
 * z) with the solute S:
 */
static real ljc (const Site *A, const Site S[], int nsites)
{
    /* Sum force field contribution from all solute sites: */
    real field = 0.0;

    for (int site = 0; site < nsites; site++) {

        /* Interaction parameters for a pair of LJ sites: */
        real e2 = sqrt (A->epsilon * S[site].epsilon);
        real s2 = 0.5 * (A->sigma + S[site].sigma);

        /* Distance from a grid point to this site: */
        real r_s = sqrt (SQR(A->x[0] - S[site].x[0]) +
                         SQR(A->x[1] - S[site].x[1]) +
                         SQR(A->x[2] - S[site].x[2]));

        /* 1. Lennard-Jones */
        field += Lennard_Jones (r_s, e2, s2);

        /* 2. Coulomb,  short range part.  For  historical reasons the
           overall scaling factor, the  product of solvent- and solute
           site charges, is handled by the function itself: */
        field += Coulomb_short (r_s, A->charge * S[site].charge);
    }

    return field;
}

/*
 * Charge  density  of  the solute  S  at  (x,  y, z).   Solvent  site
 * parameters (epsilon, sigma, charge) unused only the location of the
 * solvent site A->x is used here.
 *
 * Each gaussian is evaluated as:
 *
 *   rho(r) = q * [ G / sqrt(pi)]^3 * exp[-G^2 * (r - x0)^2]
 */
static real rho (const Site *A, const Site S[], int nsites)
{
    /* G  is predefind  in bgy3d_SolventParameters.h  FIXME:  make the
       gaussian width a property of  the (solute) site in the same way
       as the charge of the site. */
    real prefac = pow(G / sqrt(M_PI), 3.0);

    /* Sum Gaussian contributions from all solute sites: */
    real field = 0.0;

    for (int site = 0; site < nsites; site++) {

        /* Square of the distance from a grid point to this site: */
        real r2 = (SQR(A->x[0] - S[site].x[0]) +
                   SQR(A->x[1] - S[site].x[1]) +
                   SQR(A->x[2] - S[site].x[2]));

        /* Gaussian  distribution, note  that G  is not  a  width, but
           rather an inverse of it: */
        field += prefac * S[site].charge * exp(- G * G * r2);
    }

    return field;
}

/*
  Solve  Poisson  Equation  in  Fourier space  and  get  elestrostatic
  potential by inverse FFT.

  Vec uc is intent(out).
  Vec rho is intent(in).
  real q is the overall factor.

  Side effects (do not rely on them):

  Appears to  use BHD->g_fft and BHD->fft_scratch  as working storage.
  As a  matter of  fact, it  appears that one  could provide  the same
  factual parameter  for rho and  uc to effectively solve  the Poisson
  equation "in place".
*/
void poisson (BGY3dH2OData BHD, Vec uc, Vec rho, real q)
{
    int x[3], n[3], i[3], ic[3], N[3], index;
    real h[3], interval[2], k2, fac, L, h3;
    fftw_complex *fft_work;

    interval[0] = BHD->PD->interval[0];
    interval[1] = BHD->PD->interval[1];
    L = interval[1] - interval[0];
    FOR_DIM
        h[dim] = BHD->PD->h[dim];
    FOR_DIM
        N[dim] = BHD->PD->N[dim];
    h3 = h[0] * h[1] * h[2];

    /* FIXME: are we using this  as a scratch storage? Do we overwrite
       something important? */
    fft_work = BHD->g_fft;

    /* Get local portion of the grid */
    DAGetCorners(BHD->da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

    /* Get FFT of rho: rho(i, j, k) -> fft_rho(kx, ky, kz) placed into
       fft_work(kx, ky, kz): */
    ComputeFFTfromVec_fftw(BHD->da, BHD->fft_plan_fw, rho, fft_work, BHD->fft_scratch, x, n, 0);

    /*
      Solving Poisson Equation (SI units) with FFT and IFFT:

          - LAPLACIAN U (x, y, z) = (1 / epsilon0) rho(x, y, z)
                       c

      because of x = i h, y = j  h, and z = k h, with grid spacing h =
      L / n:

             2   2
          - n / L  LAPLACIAN uc(i, j, k) = (1 / epsilon0) rho(i, j, k)

      FFT (see FFTW manual "What FFTW Really Computes"):

                                    2          2    2
      fft_uc(kx, ky, kz) = 1 / [4 pi epsilon0 k  / L ] fft_rho(kx, ky, kz)

      with

           2    2    2    2
          k = kx + ky + kz

      IFFT (see FFTW manual "What FFTW Really Computes"):

      because: IFFT(fft_uc(kx, ky, kz)) = n^3 * uc(i, j, k)

                     3   3
      uc(i, j, k) = h / L  * IFFT(fft_uc(kx, ky, kz))
    */

    // EPSILON0INV = 1 / 4 * pi * epsilon0:
    real scale = q * EPSILON0INV / M_PI * h3 / (L * L * L);

    index = 0;
    for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
        for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
            for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++) {

                FOR_DIM {
                    if( i[dim] <= N[dim] / 2)
                        ic[dim] = i[dim];
                    else
                        ic[dim] = i[dim] - N[dim];
                }

                if (ic[0] == 0 && ic[1] == 0 && ic[2] == 0) {
                    /* No point to scale zeros, obviousely: */
                    fft_work[index].re = 0;
                    fft_work[index].im = 0;
                }
                else {
                    k2 = (SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0])) / SQR(L);

                    fac = scale / k2;

                    /* Here  we compute  in place:  uc(kx, ky,  kz) :=
                       scale * rho(kx, ky, kz) / k^2 */
                    fft_work[index].re = fac * fft_work[index].re;
                    fft_work[index].im = fac * fft_work[index].im;
                }
                index++;
            }
    // NOT NEEDED: VecSet(uc, 00.0);

    /* uc := IFFT(uc(kx, ky, kz)) */
    ComputeVecfromFFT_fftw(BHD->da, BHD->fft_plan_bw, uc, fft_work, BHD->fft_scratch, x, n, 0);
}
