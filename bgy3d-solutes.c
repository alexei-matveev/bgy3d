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
 * (short  range)   force  field  and  long  range   Coulomb  for  the
 * (currently) two solvent sites.
 */

#include "bgy3d.h"
#include "bgy3d-solvents.h"
#include "bgy3dH2O.h"           /* Coulomb_short() */
#include "bgy3d-fft.h"          /* ComputeFFTfromVec_fftw(),
                                   ComputeVecfromFFT_fftw() */
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

static void poisson (State *BHD, Vec uc, Vec rho, real q);

/*
 * These two functions  obey the same interface. They  are supposed to
 * get (1)  parameters of  the solvent site  such as its  location and
 * force field  parameters, and  (2) a description  of the  solute and
 * return a  real number such as  an interaction energy  or the charge
 * density:
 */
static real ljc (const Site *A, int n, const Site S[n]);
static real rho (const Site *A, int n, const Site S[n]);

/*
 * This function expects a callback obeying the above interface as one
 * of the arguments:
 */
static void field (DA da, const ProblemData *PD,
                   Site A, int n, const Site S[n],
                   real fact,
                   real (*f)(const Site *A, int n, const Site S[n]),
                   Vec v);
static void read_charge_density (DA da, const ProblemData *PD,
				const char *filename, real fact, Vec v);

// FIXME: maybe #include "solutes.h" instead?

/*********************************/
/* HCl */
/*********************************/

static const Solute HydrogenChloride =
    {"Hydrogen chloride", 2,
     {{"H", {0.6285, 0.0, 0.0}, 2.735, 0.03971, 0.2},
      {"Cl", {-0.6285, 0.0, 0.0}, 3.353, 0.51434, -0.2}}};

/*********************************/
/* CS2 */
/*********************************/

static const Solute CarbonDisulfide =
    {"Carbon disulfide", 3,
     {{"C", {0.0, 0.0, 0.0}, 3.2, 0.10128, -0.308},
      {"S1", {-1.56, 0.0, 0.0}, 3.52, 0.395, 0.154},
      {"S2", {1.56, 0.0, 0.0}, 3.52, 0.395, 0.154}}};

/*********************************/
/* Water */
/*********************************/

static const Solute Water =
    {"Water", 3,
     {{"O", {-0.2929, 0.0, 0.0}, 3.1506, 0.1521, -0.834},
      {"OH", {0.2929, 0.757, 0.0}, 0.4, 0.046, 0.417},
      {"OH", {0.2929, -0.757, 0.0}, 0.4, 0.046, 0.417}}};

/*********************************/
/* Methanol */
/*********************************/

static const Solute Methanol =
    {"Methanol", 6,
     {{"C", {-0.748, -0.015, 0.024}, 3.5, 0.066, 0.145},
      {"HC1", {-1.293, -0.202, -0.901}, 2.5, 0.03, 0.04},
      {"HC2", {-1.263, 0.754, 0.6}, 2.5, 0.03, 0.04},
      {"HC3", {-0.699, -0.934, 0.609}, 2.5, 0.03, 0.04},
      {"O", {0.558, 0.42, -0.278}, 3.12, 0.17, -0.683},
      {"OH", {0.716, 1.404, 0.137}, 0.4, 0.04, 0.418}}};

/* BUTANOIC ACID */
/* H1 sigma and epsilon adopted */

static const Solute ButanoicAcid =
    {"Butanoic Acid", 14,
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
    {"Hexane", 20,
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

/*
 * These are  the two  solvent sites.  Coordinates  will not  be used.
 * Respective parameters are #defined  elsewhere. Also do not take the
 * names  of the  sites  literally.   The same  structure  is used  to
 * represent all (2-site) solvents, such as HCl.
 *
 * FIXME: in  bgy3d_solute_field() function  below it is  assumed that
 * the number of solvent sites is exactly two:
 */
static const Site solvent[] =
  {{"h", {0.0, 0.0, 0.0}, sH, eH, qH}, /* dont use sH, eH, qH below */
   {"o", {0.0, 0.0, 0.0}, sO, eO, qO}}; /* same for sO, eO, qO */

/* Get solute sites and name by index. Functions that do the real work
   operate on array of sites: */
void bgy3d_solute_get (int solute, int *n, const Site **sites, const char **name)
{
    assert (solute >= 0 && solute <= 5);

    *n = solutes[solute]->n;
    *sites = solutes[solute]->sites;
    *name = solutes[solute]->name;
}

/*
 * Create initial solute data.
 *
 * XXX: See (5.106) and (5.08) in the thesis. Return BHD->g_ini[0] and
 *      BHD->g_ini[1], for  H and  O in that  order, (beta *  (VM_LJ +
 *      VM_coulomb_short))   and   BHD->uc[0],   BHD->uc[1]  (beta   *
 *      VM_coulomb_long), but is beta missing here?
 */

void bgy3d_solute_field (State *BHD, int n, const Site S[n], real damp, real damp_LJ)
{
  PetscPrintf(PETSC_COMM_WORLD,"Recomputing solute data with damping factor %f (damp_LJ=%f)\n", damp, damp_LJ);


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
   * Fill the force field interaction of H and O (or other) sites with
   * the solute into the respective arrays.
   *
   * We  supply ljc()  as  a  callback function  that  is supposed  to
   * compute the  interaction of  a charged LJ  solvent site  with the
   * solute.
   */
#ifndef QM
  for (int i = 0; i < 2; i++)
      field (BHD->da, BHD->PD,
             solvent[i], n, S, factor, ljc,
             BHD->g_ini[i]);
#else
  /* At  this  place the  (short  range)  Coulomb  interaction of  the
    solvent  site   with  the  solute  was   deliberately  omitted  by
    specifying zero charge of the solvent site. This effectively makes
    a point  charge (Coulomb  short + Coulomb  long) to  a distributed
    Gaussian (Coulomb long only). */

  for (int i = 0; i < 2; i++) {

    Site neutral = solvent[i];  /* dont modify the global variable */
    neutral.charge = 0.0;       /* modify a copy */

    field (BHD->da, BHD->PD,
           neutral, n, S, factor, ljc,
           BHD->g_ini[i]);
  }
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

  /* electron density file */
  size_t MAX_LEN = 260;
  char filename[MAX_LEN];

  if (bgy3d_getopt_string("--load-charge", filename, MAX_LEN)){
      read_charge_density(BHD->da, BHD->PD, filename, 1.0, v);
  }
  else {
    /* 1. Put the solute density into Vec v. Due to the inter */
    field (BHD->da, BHD->PD, point, n, S, 1.0, rho, v);
  }

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
  for (int i = 0; i < 2; i++) {
      VecSet (BHD->uc[i], 0.0);
      VecAXPY (BHD->uc[i], solvent[i].charge, v);
  }

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
                vec[k][j][i] = fact * f (&A, n, S);
            }
        }
    }
    DAVecRestoreArray (da, v, &vec);
}

/*
 * Interaction of a charged LJ site (epsilon, sigma, charge) at (x, y,
 * z) with the solute S:
 */
static real ljc (const Site *A, int n, const Site S[n])
{
    /* Sum force field contribution from all solute sites: */
    real field = 0.0;

    for (int site = 0; site < n; site++) {

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
static real rho (const Site *A, int n, const Site S[n])
{
    /* G  is predefind  in bgy3d-solvents.h  FIXME:  make the
       gaussian width a property of  the (solute) site in the same way
       as the charge of the site. */
    real prefac = pow(G / sqrt(M_PI), 3.0);

    /* Sum Gaussian contributions from all solute sites: */
    real field = 0.0;

    for (int site = 0; site < n; site++) {

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
void poisson (State *BHD, Vec uc, Vec rho, real q)
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
    ComputeFFTfromVec_fftw(BHD->da, BHD->fft_plan_fw, rho, fft_work, BHD->fft_scratch);

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
    ComputeVecfromFFT_fftw(BHD->da, BHD->fft_plan_bw, uc, fft_work, BHD->fft_scratch);
}

/* Read density file generated by GPAW calculation
 * need Bohr to scale the length values back */
#define Bohr 0.52917725750691647
static void read_charge_density (DA da, const ProblemData *PD,
				const char *filename, real fact, Vec v)
{
    PetscScalar ***vec;
    char line_buffer[BUFSIZ];
    int AtomNum; /* atom numbers in the molecule */
    real corner[3], dx[3], dy[3], dz[3];
    int GridNum[3]; /* Grid numers in each direction: [0]:x, [1]:y, [2]:z */
    real h[3];
    FILE *fp;
    int i0, j0, k0;
    int ni, nj, nk;


    fp = fopen(filename, "r");
    if (fp == NULL) {
	PetscPrintf(PETSC_COMM_WORLD, "Can not open file %s. \n", filename);
	exit(1);
    }

    PetscPrintf(PETSC_COMM_WORLD, "Reading data from %s. \n", filename);
    /* Skip the first two lines */
    for (int i = 0; i < 2; i++) {
	fgets(line_buffer, sizeof(line_buffer), fp);
    }

    /* Atom numers and corner shifts */
    fscanf(fp, "%d %lf %lf %lf", &AtomNum, &corner[0], &corner[1], &corner[2]);

    /* Grid numbers and spaces */
    FOR_DIM {
	fscanf(fp, "%d %lf %lf %lf", &GridNum[dim], &dx[dim], &dy[dim], &dz[dim]);
    }

    /* Scale the grid space */
    dx[0] *= Bohr;
    dy[1] *= Bohr;
    dz[2] *= Bohr;

    /* Allocate memory */
    int electron[AtomNum];      /* electrons of each atom */
    real zero[AtomNum]; /* = 0.0 as in ase.io.cube.write_cube, don't
                           know the meaning */
    real x[AtomNum], y[AtomNum], z[AtomNum];

    for (int i = 0; i < AtomNum; i++){
	fscanf(fp, "%d %lf %lf %lf %lf", &electron[i], &zero[i], &x[i], &y[i], &z[i]);
	x[i] *= Bohr;
	y[i] *= Bohr;
	z[i] *= Bohr;
    }

    DAGetCorners(da, &i0, &j0, &k0, &ni, &nj, &nk);

    FOR_DIM {
	h[dim] = PD->h[dim];
    }


    /* FIXME: will need interpolation if grid not match */
    if ( GridNum[0] * GridNum[1] * GridNum[2] != ni * nj * nk) {
	PetscPrintf(PETSC_COMM_WORLD, "Grid size not match!\n");
	exit(1);
    }
    else if ( fabs(dx[0] - h[0]) >= 0.001 || fabs(dy[1] - h[1]) >= 0.001 || fabs(dz[2] - h[2]) >= 0.001) {
	PetscPrintf(PETSC_COMM_WORLD, "Grid space not match!\n");
	exit(1);
    }

    DAVecGetArray(da, v, &vec);

    /* electron density scaled by Bohr^3 in python script */
    real invB3 = 1. / Bohr / Bohr / Bohr;
    for (int i = i0; i < i0 + ni; i++) {
	for (int j = j0; j < j0 + nj; j++) {
	    for (int k = k0; k < k0 + nk; k++) {
		fscanf(fp, "%lf", &vec[i][j][k]);
		vec[i][j][k] *= fact * invB3;
	    }
	}
    }

    fclose(fp);
    DAVecRestoreArray(da, v, &vec);

}

