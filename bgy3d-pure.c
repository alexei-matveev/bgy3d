/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O.c,v 1.42 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* needs Site */
#include "bgy3d-getopt.h"
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-impure.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-poisson.h"      /* laplace staff */
#include "bgy3d-vec.h"          /* bgy3d_vec_map*() */
#include "bgy3d-pure.h"
#include <complex.h>            /* after fftw.h */

static const real NORM_REG = 1.0e-1;
static const real NORM_REG2 = 1.0e-2;

/* FIXME: bgy3d-solvents.h pollutes the namespace: */
#undef sH
#undef eH
#undef qH
#undef sO
#undef eO
#undef qO

/* Allocates g[m][m] with g[j][i] being aliased to g[i][j]: */
static void create_g2 (const DA da, int m, Vec g[m][m])
{
  /* Create global vectors */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (da, &g[i][j]);
        g[j][i] = g[i][j];
      }
}

static void destroy_g2 (int m, Vec g[m][m])
{
  /* Destroy global vectors */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecDestroy (g[i][j]);
}

static State* initialize_state (const ProblemData *PD)
{
  State *BHD = (State*) malloc (sizeof (*BHD));

  BHD->PD = PD;

  PetscPrintf (PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  PetscPrintf (PETSC_COMM_WORLD, "Regularization of normalization: NORM_REG = %e\n", NORM_REG);
  PetscPrintf (PETSC_COMM_WORLD, "                                 NORM_REG2 = %e\n", NORM_REG2);
  PetscPrintf (PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf (PETSC_COMM_WORLD, "beta = %f\n", PD->beta);

  BHD->rhos[0] = PD->rho; //2.*PD->rho;
  BHD->rhos[1] = PD->rho;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &BHD->fft_mat, &BHD->da, &BHD->dc);

  /* Create global vectors. FIXME:  u2, u2_fft probably differ only by
     factors: */
  create_g2 (BHD->da, 2, BHD->u2);
  create_g2 (BHD->dc, 2, BHD->u2_fft); /* complex */
  create_g2 (BHD->da, 2, BHD->c2);
  create_g2 (BHD->da, 2, BHD->u_ini);

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      FOR_DIM
        {
          DACreateGlobalVector (BHD->da, &BHD->F[i][j][dim]);
          BHD->F[j][i][dim] = BHD->F[i][j][dim];

          DACreateGlobalVector (BHD->da, &BHD->F_l[i][j][dim]);
          BHD->F_l[j][i][dim] = BHD->F_l[i][j][dim];
        }

  FOR_DIM
    DACreateGlobalVector (BHD->da, &BHD->v[dim]);

  /* Allocate memory for fft */
  FOR_DIM
    DACreateGlobalVector (BHD->dc, &BHD->fg2_fft[dim]); /* complex */

  /* Complex scratch vector. FIXME: is it used in pure code? */
  DACreateGlobalVector (BHD->dc, &BHD->fft_scratch); /* complex */
  DACreateGlobalVector (BHD->dc, &BHD->gfg2_fft);    /* complex */

  return BHD;
}



static void finalize_state (State *BHD)
{
  MPI_Barrier( PETSC_COMM_WORLD);

  destroy_g2 (2, BHD->u2);
  destroy_g2 (2, BHD->u2_fft);
  destroy_g2 (2, BHD->u_ini);
  destroy_g2 (2, BHD->c2);

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      FOR_DIM
        {
          VecDestroy (BHD->F[i][j][dim]);
          VecDestroy (BHD->F_l[i][j][dim]);
        }

  FOR_DIM
    {
      VecDestroy (BHD->v[dim]);
      VecDestroy (BHD->fg2_fft[dim]);
    }
  VecDestroy (BHD->gfg2_fft);

  VecDestroy (BHD->fft_scratch);

#ifdef L_BOUNDARY
  MatDestroy (BHD->M);
  KSPDestroy (BHD->ksp);
#endif

  DADestroy (BHD->da);
  DADestroy (BHD->dc);
  MatDestroy (BHD->fft_mat);

  free (BHD);
}

static real LJ_repulsive (real r, real epsilon, real sigma)
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
real Coulomb_short (real r, real q2)
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

real Coulomb_short_grad (real r, real rx, real q2)
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

/* Coulomb_long  (r, 1.0) =  EPSILON0INV *  erf (G  * r)  / r  in most
   general  case.   An  extra  argument  q2  is  the  overall  scaling
   factor. */
real Coulomb_long (real r, real q2)
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

real Coulomb_long_grad (real r, real rx, real q2)
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

real Coulomb (real r, real q2)
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

real Coulomb_grad (real r, real rx, real q2)
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

/* g := exp[-(u0 + du)], with a sanity check: */
void ComputeH2O_g (Vec g, Vec u0, Vec du)
{
  real pure f (real x, real y)
  {
    const real z = exp (- (x + y));
    assert (!isinf (z));
    assert (!isnan (z));
    return z;
  }
  bgy3d_vec_map2 (g, f, u0, du);
}


/*
  Long range  pair potential Vec uc  is intent(out) here,  same as its
  FFT  transform uc_fft.  Vec  fc[] is  filled with  the corresponding
  force derived by means of FFT from the potential uc.

  No side effects.
*/
void ComputeFFTfromCoulomb (State *BHD,
                            Vec uc, Vec fc[3], /* intent(out) */
                            Vec uc_fft,    /* complex, intent(out) */
                            Vec fc_fft[3], /* complex, intent(out) */
                            real factor)
{
  int x[3], n[3], i[3];

  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real L = PD->interval[1] - PD->interval[0];

  /* Get local portion of the grid */
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

          FOR_DIM
            if (i[dim] <= N[dim] / 2)
              ic[dim] = i[dim];
            else
              ic[dim] = i[dim] - N[dim];

          if (ic[0] == 0 && ic[1] == 0 && ic[2] == 0)
            {
              uc_fft_[i[2]][i[1]][i[0]] = 0.0; /* complex */
              FOR_DIM
                fc_fft_[dim][i[2]][i[1]][i[0]] = 0; /* complex */
            }
          else
            {
              const real k2 = (SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0])) / SQR(L);
              const real fac = EPSILON0INV / M_PI / k2;
              const int sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);

              /* Potential, complex with zero imaginary part: */
              uc_fft_[i[2]][i[1]][i[0]] = factor * sign * fac * exp(- k2 * SQR(M_PI) / SQR(G));

              /* Force, imaginary: */
              FOR_DIM
                fc_fft_[dim][i[2]][i[1]][i[0]] = 2.0 * M_PI * ic[dim] / L * (I * uc_fft_[i[2]][i[1]][i[0]]);
            }
        }
  DAVecRestoreArray (BHD->dc, uc_fft, &uc_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fc_fft[dim], &fc_fft_[dim]);

  /* FFT potential */
  MatMultTranspose (BHD->fft_mat, uc_fft, uc);
  VecScale (uc, 1./L/L/L);

  /* FFT force, seems to be always requested: */
  assert (fc != NULL);
  FOR_DIM
    {
      MatMultTranspose (BHD->fft_mat, fc_fft[dim], fc[dim]);
      VecScale (fc[dim], 1./L/L/L);
    }
}


/* Precompute forces and more for a pair: */
static void pair (State *BHD,
                  const real LJ_params[3],
                  Vec f_short[3], Vec f_long[3],
                  Vec u_ini, Vec c2,
                  Vec u2, Vec u2_fft,
                  real damp, real damp_LJ)
{
  DA da;
  PetscScalar ***u_ini_;
  PetscScalar ***(f_short_[3]);
  PetscScalar ***c2_;
  real r[3], r_s, h[3], interval[2], beta, L;
  int x[3], n[3], i[3];

  /* LJ parameters of pair interaction and charge product: */
  const real epsilon = LJ_params[0];
  const real sigma = LJ_params[1];
  const real q2 = LJ_params[2];

  const ProblemData *PD = BHD->PD;
  da = BHD->da;

  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  L = PD->interval[1]-PD->interval[0];
  beta = PD->beta;

  /* FIXME: only periodic[0] == {0.0, 0.0, 0.0} appears to be used: */
  real periodic[27][3] = \
    {{0, 0, 0},
     {L, 0, 0}, {-L, 0, 0}, {0, L, 0}, {0, -L, 0}, {0, 0, L}, {0, 0, -L},
     {L, L, 0}, {-L, L, 0}, {L, -L, 0}, {-L, -L, 0},
     {L, L, L}, {-L, L, L}, {L, -L, L}, {-L, -L, L},
     {L, L, -L}, {-L, L, -L}, {L, -L, -L}, {-L, -L, -L},
     {0, L, L}, {0, -L, L}, {0, L, -L}, {0, -L, -L},
     {L, 0, L}, {-L, 0, L}, {L, 0, -L}, {-L, 0, -L}};

  /* Get local portion of the grid */
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);


   FOR_DIM
    {
      VecSet (f_short[dim], 0.0);
      VecSet (f_long[dim], 0.0);
    }
  VecSet (u_ini, 0.0);

  /*********************************************/
  /* Compute fft from Coulomb potential (long) */
  /********************************************/

  /* FFT image of Coulomb long  force in BHD->fg2_fft[3] appears to be
     ignored since overwritten. So here these are work arrays: */
  ComputeFFTfromCoulomb (BHD, u2, f_long, u2_fft, BHD->fg2_fft, q2 * damp);

  FOR_DIM
    VecAXPY (f_short[dim], 1.0, f_long[dim]);

  DAVecGetArray (da, u_ini, &u_ini_);
  DAVecGetArray (da, c2, &c2_);
  FOR_DIM
    DAVecGetArray (da, f_short[dim], &f_short_[dim]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        for (int cell = 0; cell < 1; cell++) /* FIXME: one of 27 unit cells */
          {
            /* set force vectors */
            FOR_DIM
              r[dim] = i[dim] * h[dim] + interval[0] + periodic[cell][dim];

            r_s = sqrt (SQR (r[0]) + SQR (r[1]) + SQR (r[2]));

            /* Lennard-Jones */
            u_ini_[i[2]][i[1]][i[0]] +=
              damp_LJ * beta * Lennard_Jones (r_s, epsilon, sigma);

            /* Coulomb short */
            u_ini_[i[2]][i[1]][i[0]] +=
              damp * beta * Coulomb_short (r_s, q2);

            FOR_DIM
              {
                /* Lennard-Jones */
                f_short_[dim][i[2]][i[1]][i[0]] +=
                  damp_LJ * Lennard_Jones_grad (r_s, r[dim], epsilon, sigma);

                /* Coulomb short */
                f_short_[dim][i[2]][i[1]][i[0]] +=
                  damp * Coulomb_short_grad (r_s, r[dim], q2);
              }

            /* deterministic correction */
            c2_[i[2]][i[1]][i[0]] =
              exp (-beta * LJ_repulsive (r_s, epsilon, sigma));
          }

  DAVecRestoreArray (da, u_ini, &u_ini_);
  DAVecRestoreArray (da, c2, &c2_);
  FOR_DIM
    DAVecRestoreArray (da, f_short[dim], &f_short_[dim]);
}

static void RecomputeInitialData (State *BHD, real damp, real damp_LJ)
{
  PetscPrintf (PETSC_COMM_WORLD,
               "Recomputing initial data with damping factor %f (damp_LJ=%f)\n",
               damp, damp_LJ);

  int m;                        /* number of solvent sites */
  const Site *solvent;          /* solvent[m] */

  /* Get the number of solvent sites and their parameters: */
  bgy3d_solvent_get (&m, &solvent);

  /* Original code used to print solvent params: */
  bgy3d_sites_show ("Solvent", m, solvent);

  assert (m == 2);          /* FIXME: State was allocated for m = 2 */

  /* Over all pairs: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* Pair interaction parameters: */
        real ff_params[3];
        ff_params[0] = sqrt (solvent[i].epsilon * solvent[j].epsilon);
        ff_params[1] = 0.5 * (solvent[i].sigma + solvent[j].sigma);
        ff_params[2] = solvent[i].charge * solvent[j].charge;

        pair (BHD, ff_params,
              BHD->F[i][j], BHD->F_l[i][j],
              BHD->u_ini[i][j], BHD->c2[i][j],
              BHD->u2[i][j], BHD->u2_fft[i][j],
              damp, damp_LJ);
      }
}

/* All vectors are complex here. No side effects. */
static void kapply (const State *BHD,
                    Vec fg2_fft[3],          /* intent(in) */
                    Vec g_fft, Vec coul_fft, /* intent(in) */
                    Vec dg_fft)              /* intent(out) */
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];
  const real fac = L / (2. * M_PI); /* BHD->f ist nur grad U, nicht F=-grad U  */

  complex ***g_fft_, ***dg_fft_, ***coul_fft_, ***fg2_fft_[3];
  DAVecGetArray (BHD->dc, g_fft, &g_fft_);
  DAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecGetArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* Get local portion of the grid */
  int x[3], n[3], i[3];
  DAGetCorners(BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          FOR_DIM
            if (i[dim] <= N[dim] / 2)
              ic[dim] = i[dim];
            else
              ic[dim] = i[dim] - N[dim];

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);

          complex div = 0.0;    /* complex */
          if (likely (k2 != 0))
            {
              FOR_DIM
                div += ic[dim] * fg2_fft_[dim][i[2]][i[1]][i[0]];

              const real k_fac = (h3 * fac) / k2;

              /* "I" is an imaginary unit here: */
              div *= -I * k_fac;

              /* Long  range  Coulomb part.  Note  there  is not  "-I"
                 factor here: */
              div += coul_fft_[i[2]][i[1]][i[0]];

              /* phase shift factor for x = x + L/2 */
              const int sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);

              div *=  (h3 * sign) * g_fft_[i[2]][i[1]][i[0]];
            }
          dg_fft_[i[2]][i[1]][i[0]] = div;
        }
  DAVecRestoreArray (BHD->dc, g_fft, &g_fft_);
  DAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecRestoreArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);
}

/*
  Side effects:  uses BHD->{v[], fft_scratch,  fg2_fft[], gfg2_fft} as
  work Vecs. Does (4 + 1) FFTs. One inverse.
*/
static void Compute_dg_inter (State *BHD,
                              Vec fab_s[3], Vec fab_l[3], Vec gab,
                              Vec gb,
                              Vec cab_fft, /* long range coulomb */
                              real rhob,
                              Vec dua) /* intent(out) */
{
  const ProblemData *PD = BHD->PD;
  const real L = PD->interval[1] - PD->interval[0];

  Vec gb_fft = BHD->fft_scratch;
  Vec *fg2_fft = BHD->fg2_fft;  /* fg2_fft[3] */
  Vec dua_fft = BHD->gfg2_fft;

  /************************************************/
  /* rhob * FS*gab gb */
  /************************************************/

  /*
    This computes the k-space representation of the weighted force fft
    (f_ab  g_ab).  This  force  is the  precursor  to the  convolution
    kernel K_ab. Again the long-range part is treated separately.
  */
  FOR_DIM
    {
      VecPointwiseMult (BHD->v[dim], gab, fab_s[dim]);

      /* special treatment: Coulomb long */
      VecAXPY (BHD->v[dim], -1.0, fab_l[dim]);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }

  /* The  convolution will  be  computed in  momentum  space, so  also
     compute fft (g_b): */
  MatMult (BHD->fft_mat, gb, gb_fft);

  /* This  effectively applies  K_ab  to g_b  in  k-space and  returns
     k-space du_a: */
  kapply (BHD, fg2_fft, gb_fft, cab_fft, dua_fft);

  /* Transform the result of the convolution to real space du_a: */
  MatMultTranspose (BHD->fft_mat, dua_fft, dua);

  VecScale (dua, rhob * PD->beta/L/L/L);
}


/* w := x / max(y, thresh), essentially: */
static void safe_pointwise_divide (Vec w, /* intent(out) */
                                   Vec x, /* intent(in) */
                                   Vec y, /* intent(in) */
                                   real thresh)
{
  real pure f (real x, real y)
  {
    if (y < thresh)
      return x / thresh;
    else
      return x / y;
  }
  bgy3d_vec_map2 (w, f, x, y);
}


static double pure sinc (double x)
{
  if (unlikely (x == 0.0))
    return 1.0;
  else
    return sin (x) / x;
}


/*
  Output y(k)  is the  convolution of x(k)  with ω(x)  = δ(|x| -  r) /
  4πr².  This function takes momentum  space x(k) as input and returns
  momentum space y(k).  As it  appears that time to compute sinc(k) is
  measurable, this  function operates of  array of Vecs.  These arrays
  are either of length 1 or 3, in practice.

  NOTE: appears to work for in-place transformation with x == y.

  FIXME: sinc(0) should be always 1, even at war times!
*/
static void omega (const ProblemData *PD, const DA dc,
                   const int d, Vec x_fft[d], /* intent(in) */
                   const real rab,
                   Vec y_fft[d],       /* intent(out) */
                   const real sinc0)
{
  int x[3], n[3], i[3];

  const int *N = PD->N;         /* N[3] */
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];

  /* Get local portion of the grid */
  DAGetCorners (dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***x_fft_[d], ***y_fft_[d];
  for (int p = 0; p < d; p++)
    {
      DAVecGetArray (dc, x_fft[p], &x_fft_[p]);
      DAVecGetArray (dc, y_fft[p], &y_fft_[p]);
    }

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          FOR_DIM
            if (i[dim] <= N[dim] / 2)
              ic[dim] = i[dim];
            else
              ic[dim] = i[dim] - N[dim];

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);
          const real kr = (2.0 * M_PI * rab / L) * sqrt (k2);

          /* Compute ω(k): */
          real sinc_kr;
          if (unlikely (kr == 0.0))
            sinc_kr = sinc0; /* FIXME: sinc(0) = 1, but see call sites */
          else
            sinc_kr = sinc (kr);

          const real wk = h3 * sinc_kr;

          /* Set y(k) = wk * x(k): */
          for (int p = 0; p < d; p++)
            y_fft_[p][i[2]][i[1]][i[0]] = wk * x_fft_[p][i[2]][i[1]][i[0]];
        }

  for (int p = 0; p < d; p++)
    {
      DAVecRestoreArray (dc, x_fft[p], &x_fft_[p]);
      DAVecRestoreArray (dc, y_fft[p], &y_fft_[p]);
    }
}


/*
  Compute normalization function:

     c
    n  (x) = ∫ g  (y)  g  (x - y) dy
     ab         ac      bc

  where  one   of  the   convoluted  distributions  is   a  simplified
  "rigid-bond" distribution:  ω(x) = δ(|x|  - r) / 4πr².   This latter
  distribution is fully  characterized by the distance r  which is the
  only  input parameter. Another  distribution is  supplied as  a grid
  data.    NSSA   stays   for  "normalized   site-site   superposition
  approximation"  that  is   applied  to  triple  distributions.   See
  e.g. Eq. (4.106) on p. 75 of Jager thesis.

                                      c
  Thus, given g   and r   it returns n  .
               ac      bc             ab

  This function  takes momentum space g(k) as  input.  The convolution
  result is manipulated in  the real-space rep, unfortunately.  Though
  this  manipulation is plain  screening of  the too  small (negative)
  values.

  Vec work is a complex work vector.
  Vec nab is real, intent(out).

  As a matter of  fact, the Vec work and the input  Vec gac_fft may be
  aliased. Then the input will be destroyed, of course.

  FIXME: compare the code to Compute_dg_H2O_intra_ln().
 */
static void nssa_norm_intra (const State *BHD, Vec gac_fft, real rbc,
                             Vec work, Vec nab)
{
  const real L = BHD->PD->interval[1] - BHD->PD->interval[0];

  /* Set n(k)  := ω(k) *  g(k), put result  into Vec work.   Pass both
     gac_fft and work as arrays of length 1 to omega(): */
  omega (BHD->PD, BHD->dc, 1, &gac_fft, rbc, &work, 1.0);

  /* Inverse FFT, n(k) -> n(x): */
  MatMultTranspose (BHD->fft_mat, work, nab);

  VecScale (nab, 1./L/L/L);

  /* FIXME:  make  n(x) >  0,  because it  will  appear  later in  the
     denominator. Note that the same redundant precautions are made at
     other places too. */
  real pure f (real x)
  {
    const real thresh = 1.0e-8;
    if (x < thresh)
      return thresh;
    else
      return x;
  }
  /* Transfrom in-place: */
  bgy3d_vec_map1 (nab, f, nab); /* argument aliasing! */
}


/*
  Compute  the  so-called  third,  "challenging",  term  of  the  NSSA
  expansion.  Given the input consisting of

         ~ b   c
    F  , g  , n  , and ω   represented by r
     ac   ac   ab       bc                 bc

  compute the divergence term to be added to

    u
     ab

  Long-range Coulomb is again  treated specially, thus the presence of
  the extra arguments: the long-range force fac_l and cab_fft.

  See line 3 of Eq. (4.118) p. 79 of Jager thesis.
*/
static void Compute_dg_intra (State *BHD,
                              Vec fac[3], Vec fac_l[3],
                              Vec gac, Vec nab, Vec cac_fft,
                              real rbc,
                              Vec dg, Vec dg_help)
{
  int x[3], n[3], i[3];

  const ProblemData *PD = BHD->PD;
  const DA da = BHD->da;
  const int *N = PD->N;         /* N[3] */
  Vec *fg2_fft = BHD->fg2_fft;  /* fg2_fft[3] */

  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  Vec dg_fft = BHD->gfg2_fft;
  const real L = PD->interval[1] - PD->interval[0];
  const real scale = L / (2. * M_PI); /* siehe oben ... */

  /* Get local portion of the grid */
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/

  /* fft(f1*gac) */
  FOR_DIM
    {
      VecPointwiseMult (BHD->v[dim], gac, fac[dim]);
      /* special treatment: Coulomb long */
      VecAXPY (BHD->v[dim], -1.0, fac_l[dim]);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }

  /*
    Apply omega()  in-place, note argument  aliasing.  FIXME: pretends
    that sinc(0) ==  0, so the original code. Is  it possible to claim
    that  the k  = 0  component of  weighted forces  vanishes? Because
    then, formally, zero or one would not make a difference here:
  */
  omega (BHD->PD, BHD->dc, 3, fg2_fft, rbc, fg2_fft, 0.0);

  /* int(..)/nab */
  /* Laplace^-1 * divergence */

  /*******************************************/
  /*    int/nab  */
  /*******************************************/

  /* Back  transformation of coulomb  part, divide  by nab  and forward
     transfromation */
  FOR_DIM
    {
      MatMultTranspose (BHD->fft_mat, fg2_fft[dim], BHD->v[dim]);
      VecScale (BHD->v[dim], 1./L/L/L);

      /* A   safer   version   of   VecPointwiseDivide   (BHD->v[dim],
         BHD->v[dim], nab);

         v[dim] := v[dim] / nab, essentially:
      */
      safe_pointwise_divide (BHD->v[dim], BHD->v[dim], nab, NORM_REG);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }

  /*******************************************/
  complex ***dg_fft_, ***cac_fft_;
  complex ***fg2_fft_[3];
  DAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecGetArray (BHD->dc, cac_fft, &cac_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          FOR_DIM
            if (i[dim] <= N[dim] / 2)
              ic[dim] = i[dim];
            else
              ic[dim] = i[dim] - N[dim];

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR(ic[1]) + SQR(ic[0]);
          const real kr = (2.0 * M_PI  * rbc / L) * sqrt (k2);

          /* Use fg2_fft: */
          complex div = 0.0;    /* complex */
          FOR_DIM
            div += ic[dim] * fg2_fft_[dim][i[2]][i[1]][i[0]];

          if (unlikely (kr == 0.0))
            div *= 0.0;         /* 1/k2 is undefined */
          else
            div *= -I * ((h3 * scale) / k2);

          dg_fft_[i[2]][i[1]][i[0]] = div; /* complex */

          real sinc_kr;
          if (unlikely (kr == 0.0))
            sinc_kr = 0.0; /* FIXME: sinc(0) == 1, but so the original */
          else
            sinc_kr = sinc (kr);

          /* Overwrite fg2_fft: */
          FOR_DIM
            fg2_fft_[dim][i[2]][i[1]][i[0]] = ic[dim] / scale * (I * cac_fft_[i[2]][i[1]][i[0]]) * sinc_kr;
        }
  DAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  DAVecRestoreArray (BHD->dc, cac_fft, &cac_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg_help);

  VecScale (dg_help, PD->beta/L/L/L);

  /* Back transformation  of coulomb part,  divide by nab  and forward
     transfromation */
  FOR_DIM
    {
      MatMultTranspose (BHD->fft_mat, fg2_fft[dim], BHD->v[dim]);
      VecScale (BHD->v[dim], 1./L/L/L);

      /* A   safer   version   of   VecPointwiseDivide   (BHD->v[dim],
         BHD->v[dim], nab);

         v[dim] := v[dim] / nab, essentially:
      */
      safe_pointwise_divide (BHD->v[dim], BHD->v[dim], nab, NORM_REG);

      MatMult (BHD->fft_mat, BHD->v[dim], fg2_fft[dim]);
    }

  DAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          FOR_DIM
            if (i[dim] <= N[dim] / 2)
              ic[dim] = i[dim];
            else
              ic[dim] = i[dim] - N[dim];

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);

          complex div = 0.0;    /* divergence */
          FOR_DIM
            div += ic[dim] * fg2_fft_[dim][i[2]][i[1]][i[0]];

          if (unlikely (k2 == 0))
            div *= 0.0;         /* 1/k2 is undefined */
          else
            div *= -I * ((h3 * scale) / k2);

          dg_fft_[i[2]][i[1]][i[0]] = div;
        }
  DAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg);

  VecScale (dg, PD->beta/L/L/L);

  VecAXPY (dg, 1.0, dg_help);
}


/*
  Compute   intramolecular  part.    The  first   thing  it   does  is
  transforming g(x) -> g(k). The real space representation g(x) is not
  otherwise used.

  Vec gac is intent(in).
  Vec dg is intent(out).

  Side effects: uses BHD->fft_scratch as temp Vec.

  FIXME: compare the code to nssa_norm_intra().
*/
void Compute_dg_H2O_intra_ln (State *BHD, Vec gac, real rbc, Vec dg)
{
  /* g(x) -> g(k): */
  MatMult (BHD->fft_mat, gac, BHD->fft_scratch);

  /*
    FIXME:  argument aliasing!   BHD->fft_scratch is  the  input alias
    work  array.  Its  contents  is  destroyed on  return.  Vec dg  is
    intent(out) here:
  */
  nssa_norm_intra (BHD, BHD->fft_scratch, rbc, BHD->fft_scratch, dg);

  /* -ln(g) */
  real pure f (real x)
  {
    // assert (x >= 1.0e-8);       /* See nssa_norm_intra() */
    return - log (x);
  }
  /* Transfrom in-place: */
  bgy3d_vec_map1 (dg, f, dg); /* argument aliasing! */

  /* Ensure normalization condition int(u)=0 */
  /* VecSum (dg, &k); */
  /* VecShift (dg, -k/N[0]/N[1]/N[2]); */
}


/*
  Approximate auxilary function

     c        ~ (2)              c
    Γ  (x) ~  g      = g  (x) / n  (x)
     ab        ab;c     ab       ab

  used to represent three-body distribution

                   c   b
    g    ~  g    Γ    Γ
     abc     bc   ab   ac

  in  normalized site-site superposition  approximation for  the mixed
  triplet distributions.   Of three sites the  two b, and  c belong to
  the  same rigid  species so  that their  pair distribution  is fully
  characterised by  the distance  r only.  See  e.g.  Eqs.   (4.102) -
  (4.106) on p. 75 of Jager thesis.

  Vec tab is real, intent(out).

  Side effects: used BHD->fft_scratch as work array!
*/
void bgy3d_nssa_gamma_cond (const State *BHD, Vec gac_fft, real rbc, Vec gab,
                            Vec tab) /* intent(out) */
{
  /* n(x) goes into Vec tab which is is intent(out) here: */
  nssa_norm_intra (BHD, gac_fft, rbc, BHD->fft_scratch, tab);

  /*
    t(x) =  g(x) / n(x)  (or rather t(x)  = g(x) / t(x)  with argument
    aliasing) avoiding small denominators.  Some of the commented code
    used VecPointwiseDivide() instead.
  */
  safe_pointwise_divide (tab, gab, tab, NORM_REG2); /* argument aliasing */
}


/* solve */
Vec BGY3d_solve_2site (const ProblemData *PD, Vec g_ini)
{
  const int m = 2;              /* FIXME: hardwired number of sites */
  State *BHD;
  Vec dg_new, dg_new2, work;
  real a = 0.9;
  real a1 = 0.5;
  int iter, mycount=0, upwards=0, namecount=0;
  real dg_norm[m][m];
  real dg_norm_old[m][m];

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Fourier ansatz...\n");

  BHD = initialize_state (PD);

  if (r_HH > 0)
    PetscPrintf (PETSC_COMM_WORLD, "WARNING: Solvent not a 2-Site model!\n");

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Mixing parameter: */
  const real a0 = PD->lambda;

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  Vec g[m][m], dg[m][m], t[m][m]; /* real, ij = ji */
  Vec g_fft[m][m];                /* complex, ij = ji */
  create_g2 (BHD->da, m, g);
  create_g2 (BHD->dc, m, g_fft); /* complex */
  create_g2 (BHD->da, m, dg);
  create_g2 (BHD->da, m, t);

  DACreateGlobalVector (BHD->da, &dg_new);
  DACreateGlobalVector (BHD->da, &dg_new2);
  DACreateGlobalVector (BHD->da, &work);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD->da, BHD->PD, &BHD->M, &BHD->ksp);

  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem  across iterations. Though  by now  I am  not sure  if KSP
    solver  really uses  them as  initial guess  as opposed  to always
    starting iterations from zero:
  */
  Vec x_lapl[m][m];             /* real */
  create_g2 (BHD->da, m, x_lapl);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecSet (x_lapl[i][j], 0.0);
#endif

  /* FIXME: g0[][] is a misnomer. Rename to u0[][]: */
  assert (m == 2);
  Vec (*g0)[m] = BHD->u_ini;        /* FIXME: alias! */

  /* set initial guess*/
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecSet (dg[i][j], 0.0);

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O"))
    bgy3d_vec_read2 ("dg%d%d.bin", m, dg);

  VecSet(dg_new,0.0);

  for (real damp = damp_start; damp <= 1; damp += 0.1)
    {
      RecomputeInitialData (BHD, (damp > 0.0 ? damp : 0.0), 1.0);
      PetscPrintf (PETSC_COMM_WORLD, "New lambda= %f\n", a0);

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          {
            /* Vec work is used as a temporary here: */
            bgy3d_impose_laplace_boundary (BHD, g0[i][j], work, x_lapl[i][j]);

            /* g = g0 * exp(-dg) */
            ComputeH2O_g (g[i][j], g0[i][j], dg[i][j]);

            /* Not sure if 0.0 as inital value is right. */
            dg_norm_old[i][j] = 0.0;
            dg_norm_old[j][i] = 0.0;
          }

      a1 = a0;
      a = a0;

      /* FIXME:   hardwired   distance   matrix.   Zeros   are   never
         referenced: */
      assert (m == 2);
      const real r[2][2] = {{0.0, r_HO},
                            {r_HO, 0.0}};

  for (iter = 0; iter < max_iter; iter++)
    {
      if (!(iter % 10) && iter > 0)
        a = a1;
      else if (iter == 20)
        {
          /* nothing? */
        }
      else
        a = a0;

      /* Compute FFT of g[][] for all site pairs: */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          MatMult (BHD->fft_mat, g[i][j], g_fft[i][j]);

      /* For each site pair  ij (HH, HO, and OO) compute du  in g = g0
         exp(-du): */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          {
            /* Clear accumulator, terms will be added here: */
            VecSet (dg_new, 0.0);

            /*
              Sum over  the third site k.   See the first  line of Eq.
              (4.118), p.  79 of  the Jager thesis.  This is so-called
              inter-molecular contribution. Also note that most of the
              pair quantities  are symmetric (I would say  all of them
              if I werent afraid of over-generalizations):

                Δu   = βρ  Σ  div (F   g  ) * g
                  ij     S  k       ik  ik     jk

              The site  density common for all sites  is replaced here
              by site-specific, but so far equal, density rho[k]. Note
              that in the 3-site code the sum over equivalent sites is
              treated  by scaling  the  density rho[k]  appropriately:
              e.g.  factor 2 for  H-site in H2O.  Equivalent sites are
              by  definition those  whose distributions  densities are
              equal.
            */
            for (int k = 0; k < m; k++)
              {
                Compute_dg_inter (BHD,
                                  BHD->F[i][k], BHD->F_l[i][k], g[i][k],
                                  g[j][k],
                                  BHD->u2_fft[i][k], BHD->rhos[k],
                                  dg_new2);
                VecAXPY (dg_new, 1.0, dg_new2);
              }

            /* FIXME: 3-site code does not do this: */
            VecPointwiseMult (dg_new, dg_new, BHD->c2[i][j]);

            /*
              These are two sums  over k /= i and k /=  j. See lines 2
              and 3 of Eq. (4.118), p. 79 of Jager thesis:
            */
            for (int k = 0; k < m; k++)
              {
                /*
                  Line 2  of Eq.  (4.118),  p.  79 Jager  thesis. Here
                  sites k  and i belong  to the same  solvent species.
                  Intra-molecular  correlation is  fully  described by
                  the distance r[k][i] in the rigid solvent model.

                  Note that the site density rho[k] does not enter the
                  equation explicitly.   Instead the 3-site  code used
                  the literal  factor 2 when  incrementing accumulator
                  in a call to VecAXPY():
                */
                if (k != i)
                  {
                    /* Here t is intent(out): */
                    bgy3d_nssa_gamma_cond (BHD, g_fft[i][j], r[k][i], g[j][k], t[j][k]);

                    /* Compute dg_new2 term and add to the accumulator: */
                    Compute_dg_H2O_intra_ln (BHD, t[j][k], r[k][i], dg_new2);
                    VecAXPY (dg_new, 1.0, dg_new2);
                  }

                /*
                  Line  3 of Eq.   (4.118), p.   79 Jager  thesis. The
                  so-called "numerically challenging" term. Here sites
                  k  and  j  belong   to  the  same  solvent  species.
                  Intra-molecular  correlation is  fully  described by
                  the distance r[j][k] in the rigid solvent model:
                */
                if (k != j)
                  {               /* INTRA2 */
                    /*
                      We  need  two   distinct  temporaries  to  store
                      intermediates

                       ~ (2)        k
                       g      and  n  .
                        ik;j        ij

                      These two  qualify as k  != j in the  context of
                      this summation:
                    */
                    assert (t[i][k] != t[i][j]);

                    /*
                      NOTE:  without  the  next conditional  the  next
                      computation would be redundant computation for j
                      == i.   See the call  to bgy3d_nssa_gamma_cond()
                      above:
                    */
                    if (i != j)
                      bgy3d_nssa_gamma_cond (BHD, g_fft[i][j], r[j][k], g[i][k], t[i][k]);
                    nssa_norm_intra (BHD, g_fft[i][k], r[j][k], BHD->fft_scratch, t[i][j]);
                    Compute_dg_intra (BHD,
                                      BHD->F[i][k], BHD->F_l[i][k], t[i][k],
                                      t[i][j],
                                      BHD->u2_fft[i][k], r[j][k],
                                      dg_new2, work);
                    VecAXPY (dg_new, 1.0, dg_new2);
                  }
              }

            /* Long-range Coulomb: */
            VecAXPY (dg_new, PD->beta, BHD->u2[i][j]);

            /*
              Add  an  effective  Coulomb  field of  boundary  surface
              charge  to  make  the  "potential" du  at  the  boundary
              vanish:
            */
            if (iter >= 0)
              bgy3d_impose_laplace_boundary (BHD, dg_new, work, x_lapl[i][j]);

            /*
              Mix du and du_new with a fixed ratio "a":

                du' = a * du_new + (1 - a) * du
                norm = |du_new - du|

              last arg is a temp
            */
            dg_norm[i][j] = bgy3d_vec_mix (dg[i][j], dg_new, a, work);
            dg_norm[j][i] = dg_norm[i][j];
          }

      /*
        Now  that du[]  has been  computed using  g[] of  the previous
        iteration one can safely update  g[].  Compute g := exp[-(u0 +
        du)], with a sanity check:
      */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          ComputeH2O_g (g[i][j], g0[i][j], dg[i][j]);

      /* (fancy) step size control */
      assert (m == 2);
      mycount++;
      if (((iter-1)%10) &&
          (dg_norm_old[0][0] < dg_norm[0][0] ||
           dg_norm_old[1][1] < dg_norm[1][1] ||
           dg_norm_old[0][1] < dg_norm[0][1]))
        {
          upwards = 1;
        }
      else if (iter>20 && !((iter-1)%10) && upwards==0 &&
               (dg_norm_old[0][0] < dg_norm[0][0] ||
                dg_norm_old[1][1] < dg_norm[1][1] ||
                dg_norm_old[0][1] < dg_norm[0][1]))
        {
          a1 /= 2.0;
          if (a1 < a0)
            a1 = a0;
          mycount=0;
        }
      else
        upwards=0;

      if(mycount>20 && a1<=0.5)
        {
          a1*=2.;
          mycount=0;
        }
      else if(mycount>20 && a1>0.5)
        {
          a1=1.0;
          mycount=0;
        }

      for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) /* Yes, full range! */
          dg_norm_old[i][j] = dg_norm[i][j];

      PetscPrintf (PETSC_COMM_WORLD, "%03d ", iter + 1);
      PetscPrintf (PETSC_COMM_WORLD, "a=%f ", a);
      PetscPrintf (PETSC_COMM_WORLD, "HO=%e ", dg_norm[0][1]);
      PetscPrintf (PETSC_COMM_WORLD, "HH=%e ", dg_norm[0][0]);
      PetscPrintf (PETSC_COMM_WORLD, "OO=%e ", dg_norm[1][1]);
      PetscPrintf (PETSC_COMM_WORLD, "count=%3d upwards=%1d", mycount, upwards);
      PetscPrintf (PETSC_COMM_WORLD, "\n");

      if (dg_norm[0][0] <= norm_tol &&
          dg_norm[1][1] <= norm_tol &&
          dg_norm[0][1] <= norm_tol)
          break;
    }

  /* FIXME: Debug  output from every iteration  with different overall
     scale factors damp.  Remove when no more needed. */
  {
    char fmt[20];
    snprintf (fmt, sizeof fmt, "vec%%d%%d-%d.m", namecount++);
    bgy3d_vec_save_ascii2 (fmt, m, g);
  }

  /* Save dg[][] to binary files: */
  if (bgy3d_getopt_test ("--save-H2O"))
    bgy3d_vec_save2 ("dg%d%d.bin", m, dg);

  /* Save g2[][] to binary files: */
  bgy3d_vec_save2 ("g%d%d.bin", m, g);

    }

  destroy_g2 (m, g);
  destroy_g2 (m, g_fft);
  destroy_g2 (m, t);
  destroy_g2 (m, dg);
  destroy_g2 (m, x_lapl);

  VecDestroy (dg_new);
  VecDestroy (dg_new2);
  VecDestroy (work);

  finalize_state (BHD);

  return PETSC_NULL;
}


/* solve with product ansatz g=g0*dg */
Vec BGY3d_solve_3site (const ProblemData *PD, Vec g_ini)
{
  const int m = 2; /* FIXME: hardwired number of inequivalent sites */
  State *BHD;
  Vec dg_new, dg_new2, work;
  real a = 0.9;
  real a1 = 0.5;
  int iter, mycount=0, upwards=0, namecount=0;
  real dg_norm[2][2];
  real dg_norm_old[2][2];

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Fourier ansatz...\n");

  BHD = initialize_state (PD);

  if (r_HH < 0)
    {
      PetscPrintf (PETSC_COMM_WORLD, "Solvent not a 3-Site model!\n");
      exit(1);
    }

  BHD->rhos[0] = 2.*BHD->rhos[0];

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Mixing parameter: */
  const real a0 = PD->lambda;

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  Vec g[m][m], dg[m][m], t[m][m]; /* real, ij = ji */
  Vec g_fft[m][m];                /* complex, ij = ji */
  create_g2 (BHD->da, m, g);
  create_g2 (BHD->dc, m, g_fft); /* complex */
  create_g2 (BHD->da, m, dg);
  create_g2 (BHD->da, m, t);

  DACreateGlobalVector (BHD->da, &dg_new);
  DACreateGlobalVector (BHD->da, &dg_new2);
  DACreateGlobalVector (BHD->da, &work);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD->da, BHD->PD, &BHD->M, &BHD->ksp);

  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem  across iterations. Though  by now  I am  not sure  if KSP
    solver  really uses  them as  initial guess  as opposed  to always
    starting iterations from zero:
  */
  Vec x_lapl[2][2];             /* real */
  create_g2 (BHD->da, 2, x_lapl);

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      VecSet (x_lapl[i][j], 0.0);
#endif

  Vec g0[m][m];
  g0[0][0] = BHD->u_ini[0][0];
  g0[1][1] = BHD->u_ini[1][1];
  g0[0][1] = BHD->u_ini[0][1];

  /* set initial guess*/
  VecSet(dg[0][0],0);
  VecSet(dg[1][1],0);
  VecSet(dg[0][1],0);

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O"))
    bgy3d_vec_read2 ("dg%d%d.bin", 2, dg);

  VecSet(dg_new,0.0);

  for (real damp = damp_start; damp <= damp_start; damp += 0.1)
    {
      RecomputeInitialData (BHD, (damp > 0.0 ? damp : 0.0), 1.0);
      PetscPrintf (PETSC_COMM_WORLD, "New lambda= %f\n", a0);

      /* Vec work is used as a temporary here: */
      bgy3d_impose_laplace_boundary (BHD, g0[0][0], work, x_lapl[0][0]);
      bgy3d_impose_laplace_boundary (BHD, g0[1][1], work, x_lapl[1][1]);
      bgy3d_impose_laplace_boundary (BHD, g0[0][1], work, x_lapl[0][1]);

      /* g=g0*exp(-dg) */
      ComputeH2O_g (g[0][1], g0[0][1], dg[0][1]);
      ComputeH2O_g (g[0][0], g0[0][0], dg[0][0]);
      ComputeH2O_g (g[1][1], g0[1][1], dg[1][1]);

      a1 = a0;
      a = a0;

  for (iter = 0; iter < max_iter; iter++)
    {
      if (!(iter % 10) && iter > 0)
        a = a1;
      else if (iter == 20)
        {
          /* nothing? */
        }
      else
        a = a0;

      PetscPrintf(PETSC_COMM_WORLD,"iter %d: dg function norms:\t", iter+1);

      /* Compute FFT of g[][] for all site pairs: */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          MatMult (BHD->fft_mat, g[i][j], g_fft[i][j]);

      /* f=integral(g) */
      if (1)                    /* kflg was set when with -pair */
        {
          /* g_OH */
          VecSet (dg_new, 0.0);
          Compute_dg_inter (BHD,
                            BHD->F[1][1], BHD->F_l[1][1], g[1][1],
                            g[0][1],
                            BHD->u2_fft[1][1], BHD->rhos[1],
                            dg_new2);
          VecAXPY (dg_new, 1.0, dg_new2);

          Compute_dg_inter (BHD,
                            BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                            g[0][0],
                            BHD->u2_fft[0][1], BHD->rhos[0],
                            dg_new2);
          VecAXPY (dg_new, 1.0, dg_new2);

          bgy3d_nssa_gamma_cond (BHD, g_fft[0][1], r_HO, g[0][0], t[0][0]);
          Compute_dg_H2O_intra_ln (BHD, t[0][0], r_HO, dg_new2); /* t is intent(in) */
          VecAXPY(dg_new, 2.0, dg_new2);

          /* tO = gHO/int(gHO wHH) */
          bgy3d_nssa_gamma_cond (BHD, g_fft[0][1], r_HH, g[0][1], t[1][1]);
          nssa_norm_intra (BHD, g_fft[0][1], r_HH, BHD->fft_scratch, t[0][1]);
          Compute_dg_intra (BHD, BHD->F[0][1], BHD->F_l[0][1],
                            t[1][1], t[0][1],
                            BHD->u2_fft[0][1], r_HH, dg_new2, work);
          VecAXPY(dg_new, 1.0, dg_new2);

          bgy3d_nssa_gamma_cond (BHD, g_fft[0][1], r_HO, g[1][1], t[1][1]);
          nssa_norm_intra (BHD, g_fft[1][1], r_HO, BHD->fft_scratch, t[0][1]);
          Compute_dg_intra (BHD, BHD->F[1][1], BHD->F_l[1][1],
                            t[1][1], t[0][1],
                            BHD->u2_fft[1][1], r_HO, dg_new2, work);
          VecAXPY(dg_new, 1.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[0][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, work, x_lapl[0][1]);

          dg_norm[0][1] = bgy3d_vec_mix (dg[0][1], dg_new, a, work);

          PetscPrintf (PETSC_COMM_WORLD, "HO= %e  (%f)  ", dg_norm[0][1], a);

          /* g_H */
          VecSet (dg_new, 0.0);
          Compute_dg_inter (BHD,
                            BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                            g[0][1],
                            BHD->u2_fft[0][1], BHD->rhos[1],
                            dg_new2);
          VecAXPY (dg_new, 1.0, dg_new2);

          Compute_dg_inter (BHD,
                            BHD->F[0][0], BHD->F_l[0][0], g[0][0],
                            g[0][0],
                            BHD->u2_fft[0][0], BHD->rhos[0],
                            dg_new2);
          VecAXPY (dg_new, 1.0, dg_new2);

          bgy3d_nssa_gamma_cond (BHD, g_fft[0][0], r_HH, g[0][0], t[0][0]);
          Compute_dg_H2O_intra_ln (BHD, t[0][0], r_HH, dg_new2); /* t is intent(in) */
          VecAXPY(dg_new, 1.0, dg_new2);

          bgy3d_nssa_gamma_cond (BHD, g_fft[0][0], r_HO, g[0][1], t[0][1]);
          Compute_dg_H2O_intra_ln (BHD, t[0][1], r_HO, dg_new2); /* t is intent(in) */
          VecAXPY(dg_new, 1.0, dg_new2);

          /* tO = gH/int(gH wHH) */
          bgy3d_nssa_gamma_cond (BHD, g_fft[0][0], r_HH, g[0][0], t[1][1]);
          nssa_norm_intra (BHD, g_fft[0][0], r_HH, BHD->fft_scratch, t[0][0]);
          Compute_dg_intra (BHD, BHD->F[0][0], BHD->F_l[0][0],
                            t[1][1], t[0][0],
                            BHD->u2_fft[0][0], r_HH, dg_new2, work);
          VecAXPY(dg_new, 1.0, dg_new2);

          bgy3d_nssa_gamma_cond (BHD, g_fft[0][0], r_HO, g[0][1], t[0][1]);
          nssa_norm_intra (BHD, g_fft[0][1], r_HO, BHD->fft_scratch, t[0][0]);
          Compute_dg_intra (BHD, BHD->F[0][1], BHD->F_l[0][1],
                            t[0][1], t[0][0],
                            BHD->u2_fft[0][1], r_HO, dg_new2, work);
          VecAXPY(dg_new, 1.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[0][0]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, work, x_lapl[0][0]);

          dg_norm[0][0] = bgy3d_vec_mix (dg[0][0], dg_new, a, work);

          PetscPrintf (PETSC_COMM_WORLD, "H= %e  (%f)  ", dg_norm[0][0], a);

          /* g_O */
          VecSet (dg_new, 0.0);
          Compute_dg_inter (BHD,
                            BHD->F[0][1], BHD->F_l[0][1], g[0][1],
                            g[0][1],
                            BHD->u2_fft[0][1], BHD->rhos[0],
                            dg_new2);
          VecAXPY (dg_new, 1.0, dg_new2);

          Compute_dg_inter (BHD,
                            BHD->F[1][1], BHD->F_l[1][1], g[1][1],
                            g[1][1],
                            BHD->u2_fft[1][1], BHD->rhos[1],
                            dg_new2);
          VecAXPY (dg_new, 1.0, dg_new2);


          bgy3d_nssa_gamma_cond (BHD, g_fft[1][1], r_HO, g[0][1], t[0][1]);
          Compute_dg_H2O_intra_ln (BHD, t[0][1], r_HO, dg_new2); /* t is intent(in) */
          VecAXPY(dg_new, 2.0, dg_new2);

          bgy3d_nssa_gamma_cond (BHD, g_fft[1][1], r_HO, g[0][1], t[0][1]);
          nssa_norm_intra (BHD, g_fft[0][1], r_HO, BHD->fft_scratch, t[1][1]);
          Compute_dg_intra (BHD, BHD->F[0][1], BHD->F_l[0][1],
                            t[0][1], t[1][1],
                            BHD->u2_fft[0][1], r_HO, dg_new2, work);
          VecAXPY(dg_new, 2.0, dg_new2);

          VecAXPY(dg_new, PD->beta, BHD->u2[1][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, dg_new, work, x_lapl[1][1]);

          dg_norm[1][1] = bgy3d_vec_mix (dg[1][1], dg_new, a, work);

          PetscPrintf (PETSC_COMM_WORLD, "O= %e  (%f)  ", dg_norm[1][1], a);

          /* ende: */
          ComputeH2O_g (g[0][1], g0[0][1], dg[0][1]);
          ComputeH2O_g (g[0][0], g0[0][0], dg[0][0]);
          ComputeH2O_g (g[1][1], g0[1][1], dg[1][1]);
        } /* of if (1) */

      /* (fancy) step size control */
      mycount++;
      if (((iter-1)%10) &&
          (dg_norm_old[0][0] < dg_norm[0][0] ||
           dg_norm_old[1][1] < dg_norm[1][1] ||
           dg_norm_old[0][1] < dg_norm[0][1]))
        {
          upwards = 1;
        }
      else if (iter>20 && !((iter-1)%10) && upwards==0 &&
               (dg_norm_old[0][0] < dg_norm[0][0] ||
                dg_norm_old[1][1] < dg_norm[1][1] ||
                dg_norm_old[0][1] < dg_norm[0][1]))
        {
          a1 /= 2.0;
          if (a1 < a0)
            a1 = a0;
          mycount=0;
        }
      else
        upwards=0;

      if(mycount>20 && a1<=0.5)
        {
          a1*=2.;
          mycount=0;
        }
      else if(mycount>20 && a1>0.5)
        {
          a1=1.0;
          mycount=0;
        }
      PetscPrintf(PETSC_COMM_WORLD,"count= %d  upwards= %d", mycount, upwards);
      dg_norm_old[0][0] = dg_norm[0][0];
      dg_norm_old[1][1] = dg_norm[1][1];
      dg_norm_old[0][1] = dg_norm[0][1];
      /*********************************/

      if (dg_norm[0][0] <= norm_tol &&
          dg_norm[1][1] <= norm_tol &&
          dg_norm[0][1] <= norm_tol)
        break;

      PetscPrintf(PETSC_COMM_WORLD,"\n");
    }

  /* FIXME: Debug  output from every iteration  with different overall
     scale factors damp.  Remove when no more needed. */
  {
    char fmt[20];
    snprintf (fmt, sizeof fmt, "vec%%d%%d-%d.m", namecount++);
    bgy3d_vec_save_ascii2 (fmt, m, g);
  }

  /* Save dg[][] to binary files: */
  if (bgy3d_getopt_test ("--save-H2O"))
    bgy3d_vec_save2 ("dg%d%d.bin", 2, dg);

  /* Save g2[][] to binary files: */
  bgy3d_vec_save2 ("g%d%d.bin", 2, g);

    }

  destroy_g2 (m, g);
  destroy_g2 (m, g_fft);
  destroy_g2 (m, t);
  destroy_g2 (m, dg);
  destroy_g2 (m, x_lapl);

  VecDestroy (dg_new);
  VecDestroy (dg_new2);
  VecDestroy (work);

  finalize_state (BHD);

  return PETSC_NULL;
}

