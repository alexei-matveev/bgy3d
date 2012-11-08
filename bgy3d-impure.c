/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2OS.c,v 1.20 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"

#include "bgy3d-solvents.h"
#include "bgy3d-getopt.h"
#include "bgy3d-solutes.h"
#include "bgy3d-pure.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-poisson.h"      /* laplace staff */
#include "bgy3d-impure.h"

#ifndef L_BOUNDARY_MG
#include "bgy3d-multigrid.h"    /* InitializeDMMGSolver */
#endif

#include <float.h>              /* DBL_MAX */
#include <stdbool.h>            /* bool */
#include <complex.h>            /* after fftw.h */

/*
 * These are  the two  solvent sites.  Coordinates  will not  be used.
 * Respective parameters are #defined  elsewhere. Also do not take the
 * names  of the  sites  literally.   The same  structure  is used  to
 * represent all (2-site) solvents, such as HCl.
 *
 * FIXME:  at many places  it is  assumed that  the number  of solvent
 * sites is exactly two:
 */
static const Site solvent[] =
  {{"h", {0.0, 0.0, 0.0}, sH, eH, qH}, /* dont use sH, eH, qH below */
   {"o", {0.0, 0.0, 0.0}, sO, eO, qO}}; /* same for sO, eO, qO */

#undef sH
#undef eH
#undef qH
#undef sO
#undef eO
#undef qO

static void load (const State *BHD, Vec g2[2][2]);

static State initialize_state (const ProblemData *PD)
{
  State BHD;

  BHD.PD = PD;

  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", PD->interval[0], PD->interval[1]);
  //PetscPrintf(PETSC_COMM_WORLD, "Boundary smoothing parameters : SL= %f  SR= %f\n", SL, SR);
  //PetscPrintf(PETSC_COMM_WORLD, "ZEROPAD= %f\n", ZEROPAD);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", PD->h[0]);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", PD->beta);
  /******************************/
  BHD.rhos[0] = PD->rho;
  BHD.rhos[1] = PD->rho;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &BHD.fft_mat, &BHD.da, &BHD.dc);

#ifdef L_BOUNDARY_MG
  /* multigrid, apparently needs two descriptors: */
#error "Need BHD.da_dmmg"
#endif

  const DA da = BHD.da;         /* shorter alias */

  /* Create global vectors */
  DACreateGlobalVector (da, &BHD.g_ini[0]);
  DACreateGlobalVector (da, &BHD.g_ini[1]);

  BHD.gHO_ini = PETSC_NULL;     /* unused with impurities */

  DACreateGlobalVector (da, &BHD.u2[0][0]);
  DACreateGlobalVector (da, &BHD.u2[1][1]);
  DACreateGlobalVector (da, &BHD.u2[0][1]);

  /* Pair quantities  here, use symmetry wrt  (i <-> j)  to save space
     and work: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        /* FIXME: see bgy3d_load_vec below, arent we overwriting these
           g2 vecs? */
        DACreateGlobalVector (da, &BHD.g2[i][j]);
        BHD.g2[j][i] = BHD.g2[i][j];

        /* Used with pure solvent only: */
        FOR_DIM
          {
            BHD.F[i][j][dim] = PETSC_NULL;
            BHD.F[j][i][dim] = PETSC_NULL;

            BHD.F_l[i][j][dim] = PETSC_NULL;
            BHD.F_l[j][i][dim] = PETSC_NULL;
          }
      }

  BHD.pre = PETSC_NULL;         /* used for newton solver only */

  FOR_DIM
    DACreateGlobalVector (da, &BHD.v[dim]);


#ifdef L_BOUNDARY
  DACreateGlobalVector(da, &BHD.x_lapl[0]);
  DACreateGlobalVector(da, &BHD.x_lapl[1]);
  VecSet(BHD.x_lapl[0], 0.0);
  VecSet(BHD.x_lapl[1], 0.0);
#endif

  /* Allocate memory for fft */
  FOR_DIM
    {
      for (int i = 0; i < 2; i++)
        for (int j = 0; j <= i; j++)
          {
            DACreateGlobalVector (BHD.dc, &BHD.fs_g2_fft[i][j][dim]);
            BHD.fs_g2_fft[j][i][dim] = BHD.fs_g2_fft[i][j][dim];

            DACreateGlobalVector (BHD.dc, &BHD.fl_g2_fft[i][j][dim]);
            BHD.fl_g2_fft[j][i][dim] = BHD.fl_g2_fft[i][j][dim];
          }

      DACreateGlobalVector (BHD.dc, &BHD.fg2_fft[dim]); /* used by ComputeFFTfromCoulomb() */
    }

  /* Complex scratch vector: */
  DACreateGlobalVector (BHD.dc, &BHD.fft_scratch);

  /* FIXME: these probably differ only by factors: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD.dc, &BHD.u2_fft[i][j]);
        BHD.u2_fft[j][i] = BHD.u2_fft[i][j];
      }

  BHD.gfg2_fft = NULL;          /* not used with impurities */

  /* Not used with impurities: */
  BHD.LJ_paramsH[0] = -1;
  BHD.LJ_paramsH[1] = -1;
  BHD.LJ_paramsH[2] = -1;
  BHD.LJ_paramsO[0] = -1;
  BHD.LJ_paramsO[1] = -1;
  BHD.LJ_paramsO[2] = -1;
  BHD.LJ_paramsHO[0] = -1;
  BHD.LJ_paramsHO[1] = -1;
  BHD.LJ_paramsHO[2] = -1;

  /* FIXME: broken, see the comments inside the function itself: */
  load (&BHD, BHD.g2);

  return BHD;
}


/* This fills  an array of preallocated Vecs  with solven-solvent pair
   distributions:  */
static void load (const State *BHD, Vec g2[2][2])
{
#ifdef CS2
  ReadPairDistribution (BHD, "g11.txt", g2[1][1]);
  ReadPairDistribution (BHD, "g00.txt", g2[0][0]);
  ReadPairDistribution (BHD, "g01.txt", g2[0][1]);
#else
  (void) BHD;                   /* unused */
  /* Read  g^2 from  file  into a  pre-allocated global  (distributed)
   * vector.   Note in  bgy3d_load_vec(), a  new 'local'  Vec  will be
   * created  but the  new local  vec  isn't derived  from the  global
   * vector  here.  Since  global vectors  are used  in this  code, we
   * might  simply retain  operating  with global  vectors instead  of
   * converting and assembling between global and local vectors? */
  bgy3d_read_vec ("g00.bin", g2[0][0]);
  bgy3d_read_vec ("g01.bin", g2[0][1]);
  bgy3d_read_vec ("g11.bin", g2[1][1]);
#endif
  assert (g2[1][0] == g2[0][1]);
}


static void finalize_state (State *BHD)
{
  MPI_Barrier( PETSC_COMM_WORLD);

  /* Pair quantities here: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        VecDestroy (BHD->g2[i][j]);

        FOR_DIM
          {
            /* Used with pure solvent only: */
            assert (BHD->F[i][j][dim] == PETSC_NULL);
            assert (BHD->F_l[i][j][dim] == PETSC_NULL);

            VecDestroy (BHD->fs_g2_fft[i][j][dim]);
            VecDestroy (BHD->fl_g2_fft[i][j][dim]);
          }
      }

  FOR_DIM
    {
      VecDestroy (BHD->v[dim]);
      VecDestroy (BHD->fg2_fft[dim]);
    }

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      VecDestroy (BHD->u2_fft[i][j]);

  assert (BHD->gfg2_fft == NULL); /* not used with impurities */

  VecDestroy (BHD->fft_scratch);

  VecDestroy(BHD->g_ini[0]);
  VecDestroy(BHD->g_ini[1]);
  assert (BHD->gHO_ini == PETSC_NULL); /* unused for impurities */
  VecDestroy(BHD->u2[0][0]);
  VecDestroy(BHD->u2[1][1]);
  VecDestroy(BHD->u2[0][1]);
  assert (BHD->pre == PETSC_NULL); /* used for newton solver */

#ifdef L_BOUNDARY
  MatDestroy (BHD->M);
  KSPDestroy(BHD->ksp);
  VecDestroy(BHD->x_lapl[1]);
  VecDestroy(BHD->x_lapl[0]);
#endif
#ifdef L_BOUNDARY_MG
  DMMGDestroy(BHD->dmmg);
#endif

  DADestroy (BHD->da);
  DADestroy (BHD->dc);
  MatDestroy (BHD->fft_mat);
}


/* Fills Vec  g2 with  3D distribution derived  from the 1D  g(r) data
   from the disk.  Here g2 should be a valid allocated vector. */
void ReadPairDistribution (const State *BHD, const char *filename, Vec g2)
{
  DA da;
  FILE *fp;
  real *xg, *g;
  real r[3], r_s, h[3], interval[2];
  int index=0;
  int x[3], n[3], i[3], k;
  PetscScalar ***g2_vec;

  da = BHD->da;
  const ProblemData *PD = BHD->PD;
  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  /* read file */
  fp = fopen(filename, "r");
  if(fp==NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Could not open file %s.\n", filename);
      exit(1);
    }
  xg= (real*) malloc(sizeof(*xg));
  g= (real*) malloc(sizeof(*g));

  while( fscanf(fp,"%lf %lf", &xg[index], &g[index]) == 2)
    {
      index++;
      xg= (real*) realloc(xg, (index+1)*sizeof(*xg));
      g= (real*) realloc(g, (index+1)*sizeof(*g));
    }
  index--;
  PetscPrintf(PETSC_COMM_WORLD,"Read %d lines from file %s, x=[%f,%f]\n", index+1, filename,
              xg[0], xg[index]);
  fclose(fp);


  /* interpolate to 3d grid */

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  DAVecGetArray(da, g2, &g2_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          FOR_DIM
            r[dim] = i[dim]*h[dim]+interval[0];
          r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

          /* find x in array */
          for(k=0; k<=index; k++)
            if(r_s<xg[k])
              break;
          if(k==0)
            g2_vec[i[2]][i[1]][i[0]] = 0;
          else if(k>=index)
            g2_vec[i[2]][i[1]][i[0]] = 1.0;
          else
            g2_vec[i[2]][i[1]][i[0]] = g[k] + (r_s-xg[k])*(g[k+1]-g[k])/(xg[k+1]-xg[k]);
          if(g2_vec[i[2]][i[1]][i[0]]<0)
            g2_vec[i[2]][i[1]][i[0]]=0;
        }
  DAVecRestoreArray(da, g2, &g2_vec);

  free(xg);
  free(g);
}

static void  pair (State *BHD,
                   const real LJ_params[3],
                   const Vec g2,
                   Vec f_short[3], Vec f_long[3], /* work arrays */
                   Vec fs_g2_fft[3],              /* intent(out) */
                   Vec fl_g2_fft[3],              /* intent(out) */
                   Vec u2, Vec u2_fft,            /* intent(out) */
                   real damp, real damp_LJ)
{
  const real epsilon = LJ_params[0]; /* geometric average of two */
  const real sigma = LJ_params[1];   /* arithmetic average of two */
  const real q2 = LJ_params[2];      /* charge product */

  const ProblemData *PD = BHD->PD;
  const real *h = PD->h;        /* h[3] */
  const DA da = BHD->da;
  const real off = PD->interval[0];

  /*
    Compute Coulomb from fft part.

    Here Vec u2  and a complex array u2_fft[]  both are intent(out) in
    the next call.  The Vec f_long[], intent(out), is  filled with the
    corresponding force.   Performs 4 FFTs.  Again note that  the only
    difference for all u2[i][j] and their FFT transform is the overall
    scaling factor  q[i] * q[j].  FIXME: why  keeping O(m^2) versions,
    with m being number of solvent sites, of almost the same field and
    repeating unnecessary FFTs?
  */
  ComputeFFTfromCoulomb (BHD, u2, f_long, u2_fft, q2);

  /*
    Sort-range  potential/force is  specific  for each  pair, on  the
    other hand:
  */
  PetscScalar ***fs_vec[3];
  FOR_DIM
    DAVecGetArray (da, f_short[dim], &fs_vec[dim]);

  /* Get local portion of the grid */
  int x[3], n[3], i[3];
  DAGetCorners (da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          real r[3], r_s;
          /* set force vectors */
          FOR_DIM
            r[dim] = i[dim] * h[dim] + off; /* FIXME: offset */

          r_s = sqrt (SQR(r[0]) + SQR(r[1]) + SQR(r[2]));

          /* Lennard-Jones + Coulomb short */
          FOR_DIM
            fs_vec[dim][i[2]][i[1]][i[0]] =
            damp_LJ * Lennard_Jones_grad (r_s, r[dim], epsilon, sigma) +
            damp * Coulomb_short_grad (r_s, r[dim], q2);
        }

  FOR_DIM
    DAVecRestoreArray (da, f_short[dim], &fs_vec[dim]);

  /*
    Compute FFT(F * g^2).

    F * g2 = (F_LJ + F_coulomb_short) * g2
           + (F_coulomb_long * g2 - F_coulomb_long)
           + F_coulomb_long

    see (5.101) and (5.102) in Jager's thesis. FFT(F_coulomb_long) has
    been calculated  as u2_fft by  ComputeFFTfromCoulomb() above.  The
    code needs at least one work vector, use this:
  */
  Vec work = BHD->v[0];

  FOR_DIM
    {
      /* First (F_LJ + F_coulomb_short) * g2: */
      VecPointwiseMult (work, g2, f_short[dim]);

      /* Next FFT((F_LJ + F_coulomb_short) * g2): */
      MatMult (BHD->fft_mat, work, fs_g2_fft[dim]);

      /* Now Coulomb long. F_coulomb_long * g2: */
      VecPointwiseMult (work, g2, f_long[dim]);

      /* Next F_coulomb_long * g2 - F_coulomb_long: */
      VecAXPY (work, -1.0, f_long[dim]);

      /* Finally FFT(F_coulomb_long * g2 - F_coulomb_long): */
      MatMult (BHD->fft_mat, work, fl_g2_fft[dim]);
    }
}

void RecomputeInitialFFTs (State *BHD, real damp, real damp_LJ)
{
  real ff_params[3];

  PetscPrintf (PETSC_COMM_WORLD,
               "Recomputing FFT data with damping factor %f (damp_LJ=%f)\n",
               damp, damp_LJ);

  /* FIXME: maybe use BHD->v[3] for one of them? */
  Vec force_short[3];           /* work vectors for pair() */
  Vec force_long[3];            /* work vectors for pair() */
  FOR_DIM
    {
      DACreateGlobalVector (BHD->da, &force_short[dim]);
      DACreateGlobalVector (BHD->da, &force_long[dim]);
    }


  /* Over all distinct solvent site pairs: */
  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        /* For   debug  only,   check   the  symmetry   of  the   pair
           quantities: */
        if (1)
          {
            assert (BHD->g2[j][i] == BHD->g2[i][j]);
            FOR_DIM
              {
                assert (BHD->fs_g2_fft[j][i][dim] == BHD->fs_g2_fft[i][j][dim]);
                assert (BHD->fl_g2_fft[j][i][dim] == BHD->fl_g2_fft[i][j][dim]);
              }
          }

        /* Pair interaction parameters: */
        ff_params[0] = sqrt (solvent[i].epsilon * solvent[j].epsilon);
        ff_params[1] = 0.5 * (solvent[i].sigma + solvent[j].sigma);
        ff_params[2] = solvent[i].charge * solvent[j].charge;

        /* Does real work: */
        pair (BHD, ff_params,
              BHD->g2[i][j],
              force_short, force_long, /* work vectors*/
              BHD->fs_g2_fft[i][j], BHD->fl_g2_fft[i][j],
              BHD->u2[j][i], BHD->u2_fft[j][i], /* FIXME: ij ji */
              damp, damp_LJ);
      }

  /* Clean up and exit: */
  FOR_DIM
    {
      VecDestroy (force_short[dim]);
      VecDestroy (force_long[dim]);
    }
}

/*
 * This is supposed to compute the k-grid representation of the BGY3dM
 * equation:
 *
 *     Δu = - β ρ div((F g2) * g)
 *
 * with "*" being  the convolution that in Fourier  rep degenerates to
 * pointwise  product,  g2 and  g  being  the (fixed)  solvent-solvent
 * site-site  distribution  function and  (the  unknown) solvent  site
 * distribution around the solute.
 *
 * The kernel dfg(kx,  ky, kz) includes the Poisson factor  1 / k2 but
 * does not include the facotor (rho * beta):
 *
 *     u(kx, ky, kz) = β ρ dfg(kx, ky, kz) g(kx, ky, kz)
 *
 * The three fg arrays are intent(in). The dfg array is intent(out).
 *
 * No known side effects.
 */
static void kernel (const DA dc,
                    const ProblemData *PD,
                    Vec fg[3],  /* complex, intent(in) */
                    Vec coul,   /* complex, intent(in) */
                    Vec dfg)    /* complex, intent(out) */
{
  int x[3], n[3], i[3], N[3], ic[3];

  FOR_DIM
    N[dim] = PD->N[dim];

  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];
  const real L3 = L * L * L;
  const real fac = L / (2.0 * M_PI); /* BHD->f ist nur grad U, nicht F=-grad U  */
  const real scale = fac / L3;

  /* Get local portion of the grid */
  DAGetCorners (dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* Loop over local portion of grid: */
  complex ***fg_[3], ***dfg_, ***coul_;
  DAVecGetArray (dc, dfg, &dfg_);
  if (coul)
    DAVecGetArray (dc, coul, &coul_);
  FOR_DIM
    DAVecGetArray (dc, fg[dim], &fg_[dim]);

  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          FOR_DIM
            {
              if (i[dim] <= N[dim] / 2)
                ic[dim] = i[dim];
              else
                ic[dim] = i[dim] - N[dim];
            }

          /* phase shift factor for x=x+L/2 */
          real sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);

          real k2 = SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0]);

          real k_fac;
          if (k2 > 0)
              k_fac = sign * h3 * h3 * scale / k2;
          else
              k_fac = 0.0;

          /*
           * Compute the  (Fourier transform of) of  the divergence of
           * the "weighted force"  vector div (F g) which  serves as a
           * convolution  kernel  in BGY3dM  equations.  Note the  the
           * derivative in momentum space involves a factor -i so that
           * real- and imaginary  parts of divergence are proportional
           * to
           *
           *     Re = + dot(k, Im(F g))
           *     Im = - dot(k, Re(F g))
           *
           * The  additional factor  -k^(-2)  effectively included  in
           * k_fac recovers the Fourier transform of the corresponding
           * Poisson solution.
           */

          /*
            Complex  arithmetics here.   "I" is  a macro  expanding to
            imaginary unit:
          */
          complex sum = 0.0;
          for (int p = 0; p < 3; p++)
            sum += ic[p] * fg_[p][i[2]][i[1]][i[0]];

          dfg_[i[2]][i[1]][i[0]] = k_fac * (-I * sum);

          /*
           * FIXME: Origin of  this occasional addition needs some
           * explanation.
           *
           * The   difference   between   Compute_H2O_interS_C()   and
           * Compute_H2O_interS() of the original code reduces to this
           * addendum.   Supply   a   NULL   pointer   for   coul   in
           * Compute_H2O_interS_C()   to    get   the   behaviour   of
           * Compute_H2O_interS().
           *
           * Long range Coulomb part (right one):
           */
          if (coul)
            dfg_[i[2]][i[1]][i[0]] += (h3 / L3) * sign * coul_[i[2]][i[1]][i[0]];
        }
  DAVecRestoreArray (dc, dfg, &dfg_);
  if (coul)
    DAVecRestoreArray (dc, coul, &coul_);
  FOR_DIM
    DAVecRestoreArray (dc, fg[dim], &fg_[dim]);
}

/*
 * This applies the kernel compured by  kernel() to FFT of g to obtain
 * an increment to "dg". The latter probably needs a better name. Dont
 * forget to clear "dg" early enough.
 *
 * Complex Vec dg is intent(inout).
 */
static void apply (const DA dc,
                   Vec ker,          /* kernel, intent(in) */
                   Vec g,            /* current g, intent(in) */
                   const real scale, /* overall scale */
                   Vec dg)           /* incremented, intent(inout) */
{
  int x[3], n[3];

  /* Get local portion of the grid */
  DAGetCorners (dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***ker_, ***g_, ***dg_;
  DAVecGetArray (dc, ker, &ker_);
  DAVecGetArray (dc, dg, &dg_);
  DAVecGetArray (dc, g, &g_);

  /* loop over local portion of grid */
  for (int k = x[2]; k < x[2] + n[2]; k++)
    for (int j = x[1]; j < x[1] + n[1]; j++)
      for (int i = x[0]; i < x[0] + n[0]; i++)
        {
          /*
           Retrive  the   precomuted  Fourier  transform   of  of  the
           divergence of  the "weighted force" vector div  (F g) which
           serves as  a convolution  kernel in BGY3dM  equations.  The
           additional  factor  -k^(-2)  effectively  included  in  the
           kernel recovers the  Fourier transform of the corresponding
           Poisson solution.   The factor scale = βρ  is not included,
           on the other hand. See kernel() for details:
           */
          dg_[k][j][i] += scale * (ker_[k][j][i] * g_[k][j][i]); /* complex */
        }
  DAVecRestoreArray (dc, ker, &ker_);
  DAVecRestoreArray (dc, dg, &dg_);
  DAVecRestoreArray (dc, g, &g_);
}

/*
 * Side effects: none, but consider efficiency.
 */
static void Compute_H2O_interS_C (const State *BHD,
                                  Vec fg2_fft[3], /* complex, intent(in) */
                                  Vec g,        /* real, intent(in) */
                                  Vec coul_fft, /* complex, intent(in) */
                                  real rho, Vec dg)
{
  /* Avoid separate VecScale at the end: */
  const real scale = rho * BHD->PD->beta;

  /************************************************/
  /* rho*F*g^2 g*/
  /************************************************/

  /* FIXME: move allocations out of the loop: */
  Vec ker_fft, g_fft, dg_fft;
  DACreateGlobalVector (BHD->dc, &ker_fft);
  DACreateGlobalVector (BHD->dc, &g_fft);
  DACreateGlobalVector (BHD->dc, &dg_fft);

  /*
    FIXME:  Move   computation  of  the  kernel  out   of  the  BGY3dM
    iterations.   Here we put  the fft  of the  kernel into  local Vec
    ker_fft:
  */
  kernel (BHD->dc, BHD->PD, fg2_fft, coul_fft,
          ker_fft);    /* result */

  /* fft(g) */
  MatMult (BHD->fft_mat, g, g_fft);

  /* Will be incremented: */
  VecSet (dg_fft, 0.0);

  /* Apply the kernel, Put result into the complex temp Vec dg_fft: */
  apply (BHD->dc, ker_fft, g_fft, scale, dg_fft);

  /* ifft(dg) */
  MatMultTranspose (BHD->fft_mat, dg_fft, dg);

  VecDestroy (ker_fft);
  VecDestroy (g_fft);
  VecDestroy (dg_fft);
}

/*
 * Side effects: none, but see Compute_H2O_interS_C().
 */
void Compute_H2O_interS (const State *BHD, /* NOTE: modifies BHD->fft dynamic arrays */
                         Vec fg2_fft[3],   /* complex, intent(in) */
                         Vec g,            /* real, intent(in) */
                         real rho, Vec dg_help)
{
    Compute_H2O_interS_C(BHD, fg2_fft, g, NULL, rho, dg_help);
}

static real ComputeCharge (const ProblemData *PD,
                           int m, const Site solvent[m],
                           const Vec g[m],
                           Vec work)
{
  const real dV = PD->h[0] * PD->h[1] * PD->h[2];
  const real rho = PD->rho;

  real total = 0.0;
  for (int i = 0; i < m; i++)
    {
      real sum;
      VecCopy (g[i], work);
      VecShift (work, -1.0);
      VecSum (work, &sum);
      total += sum * rho * dV * solvent[i].charge;
    }

  return total;
}

static void ComputedgFromg (Vec dg, Vec g0, Vec g)
{
  int local_size, i;
  PetscScalar *dg_vec, *g_vec, *g0_vec;
  real k;

  VecGetLocalSize(dg, &local_size);

  VecGetArray(dg, &dg_vec);
  VecGetArray(g0, &g0_vec);
  VecGetArray(g, &g_vec);
  for(i=0; i<local_size; i++)
    {
      k = g_vec[i];
      if(k<1.0e-8)
        k=1.0e-8;
      dg_vec[i] = -g0_vec[i] - log(k);

    }

  VecRestoreArray(dg, &dg_vec);
  VecRestoreArray(g0, &g0_vec);
  VecRestoreArray(g, &g_vec);
}

/*
 * Does the mixing:
 *
 *     dg := a * dg_new + (1 - a) * dg
 *
 * Returns the norm of the difference |dg_new - dg|.
 */
static real mix (Vec dg, Vec dg_new, real a, Vec work)
{
    real norm;

    /* Move dg */
    VecCopy(dg, work);

    /* dg' = a * dg_new + (1 - a) * dg */
    VecAXPBY(dg, a, (1-a), dg_new);

    /* work = dg - dg' = a * (dg - dg_new)

       That is why the divison by "a" below */
    VecAXPY(work, -1.0, dg);

    /* Norm of the change: */
    VecNorm(work, NORM_INFINITY, &norm);

    /* FIXME: why not computing |dg_new - dg| directly? Would be valid
       also for a == 0: */
    return norm / a;
}

/* Returns most  negative number for  zero sized arrays.   Will return
   NaN if there is any in the array. */
static double maxval (size_t n, const double x[n])
{
  double max = -DBL_MAX;
  for (size_t i = 0; i < n && !isnan (max); i++)
    {
      /* If x[i] is  NaN comparison will fail and  max will become NaN
         too:  */
      max = x[i] < max ? max : x[i];
      /* We  need  to  break  out  otherwise  in  the  next  iteration
         comparison will  fail again and  max may be overwritten  by a
         valid number ... */
    }
  return max;
}

/*
  Function to  generate solvent field that  is supposed to  act on the
  solute  electrons. At  the moment  only the  electrostatic  field is
  computed.  This  "reaction"  field  should be  consistent  with  the
  "action" of  the solute electrons  on the solvent sites  computed in
  bgy3d-solvent.c that is also pure electrostatics:
*/
static void bgy3d_solvent_field (const State *BHD, /* intent(in) */
                                 int m, const Site solvent[m], /* m == 2*/
                                 Vec g[m], /* intent(in) */
                                 Vec ve)   /* intent(out) */
{
  /* FIXME: this  code assumes  the same density  rho for  all solvent
     particles. */

  /* Reset its value to zero */
  VecSet (ve, 0.0);

  for (int i = 0; i < m; i++)
    VecAXPY (ve, solvent[i].charge, g[i]);

  /* Solve Poisson  equation in place. Note that  the output potential
     is in kcals as all energies in this code are: */
  bgy3d_poisson (BHD, ve, ve, BHD->PD->rho); /* aliasing */

  /*
   Solving Poisson  equation by  FFT results in  a potential  with the
   mean value zero over the whole volume. That is how the ve(k) is set
   for k = 0. The following code adds a correction to that field which
   ensures  that  the  potential   is  zero  at  the  boundary.   This
   corresponds  to  a  boundary  as  a grounded  metallic  cage.   The
   corrected potential  is, in effect, a superposition  of the solvent
   electrostatic field and a surface charge on that metallic cage:
  */
  {
    Vec x, work;
    DACreateGlobalVector (BHD->da, &x);
    DACreateGlobalVector (BHD->da, &work);

    bgy3d_impose_laplace_boundary (BHD, ve, work, x);

    VecDestroy (x);
    VecDestroy (work);
  }

  bgy3d_save_vec ("ve.bin", ve); /* for debugging only */
}


/*
  This function is the main entry  point for the BGY3dM equation for a
  2-site solvent and an arbitrary solute.  The two vectors in

  Vec g[2], intent(out)

  are  initialized as global  distributed arrays  and filled  with the
  solvent site  distributions. It is the responsibility  of the caller
  to destroy them when no more needed.
 */
void bgy3d_solve_with_solute (const ProblemData *PD,
                              int n, const Site solute[n],
                              void (*density)(int k, const real x[k][3], real rho[k]),
                              Vec g[2])
{
  Vec t_vec;                 /* used for all sites */
  Vec uc;                    /* Coulomb long, common for all sites. */
  Vec du[2], du_acc, work;
  Vec ve;            /* ve for solvent electrostaic potential field */
  PetscScalar du_norm[2];
  int namecount = 0;
  char nameH[20], nameO[20];

  PetscPrintf (PETSC_COMM_WORLD, "Solving BGY3dM (2-site) equation ...\n");

  State BHD = initialize_state (PD);

  if (r_HH > 0.0)
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

  /* Inverse temperature: */
  const real beta = PD->beta;

  PetscPrintf (PETSC_COMM_WORLD, "New lambda= %f\n", a0);

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD.da, BHD.PD, &BHD.M, &BHD.ksp);
#endif
#ifdef L_BOUNDARY_MG
  InitializeDMMGSolver(&BHD);
#endif

  /*
    These complex  vectors will  hold FFT of  the current  g. Allocate
    enough to  hold a  local portion  of the grid  and free  after the
    loop.
  */
  Vec g_fft[2];

  for (int i = 0; i < 2; i++)
    {
      DACreateGlobalVector (BHD.dc, &g_fft[i]); /* complex */

      DACreateGlobalVector (BHD.da, &du[i]); /* real */

      /* Here the storage for the output is allocated, the caller will
         have to destroy them: */
      DACreateGlobalVector (BHD.da, &g[i]); /* real */
    }

  /* These are  the (four)  kernels HH, HO,  OH, OO stored  as complex
     vectors in momentum space.  Note that HO = OH. */
  Vec ker_fft_S[2][2];
  Vec ker_fft_L[2][2];

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      {
        DACreateGlobalVector (BHD.dc, &ker_fft_S[i][j]); /* complex */
        ker_fft_S[j][i] = ker_fft_S[i][j];

        DACreateGlobalVector (BHD.dc, &ker_fft_L[i][j]); /* complex */
        ker_fft_L[j][i] = ker_fft_L[i][j];
      }

  /*
    There is no  point to transform each contribution  computed in the
    momentum space back, accumulate them on the k-grid in this complex
    Vec:
  */
  Vec du_acc_fft;
  DACreateGlobalVector (BHD.dc, &du_acc_fft); /* complex */

  DACreateGlobalVector (BHD.da, &du_acc);
  DACreateGlobalVector (BHD.da, &work);
  DACreateGlobalVector (BHD.da, &t_vec); /* used for all sites */
  DACreateGlobalVector (BHD.da, &uc);    /* common for all sites */
  DACreateGlobalVector (BHD.da, &ve); /* solvent electrostatic potential */

  /*
    Later u0  = beta *  (VM_LJ + VM_coulomb_short), which  is -log(g0)
    actually.  See: (5.106)  and (5.108) in Jager's thesis.  It is not
    filled with data yet, I assume.
  */
  Vec *u0 = BHD.g_ini;          /* FIXME: aliasing! */

  /* set initial guess*/
  VecSet (du[0], 0.0);
  VecSet (du[1], 0.0);

  if (bgy3d_getopt_test ("--from-g2"))
    {
      ComputedgFromg (du[0], u0[0], BHD.g2[0][1]);
      ComputedgFromg (du[1], u0[1], BHD.g2[1][1]);
    }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O"))
    {
      PetscPrintf (PETSC_COMM_WORLD, "Loading binary files...");
      du[0] = bgy3d_load_vec ("du0.bin"); /* duH */
      du[1] = bgy3d_load_vec ("du1.bin"); /* duO */
      PetscPrintf (PETSC_COMM_WORLD, "done.\n");
    }

  for (real damp = damp_start; damp <= 1; damp += 0.1)
    {

      /*
        FIXME: I guess the logic with damping factors can be made more
        straightforward:
      */
      real damp_LJ = (damp >= 0 ? 1.0 : 0.0); /* yes, >=, so the
                                                 original */

      /*
        Return F  * g2.  Note the  calculation of F is  divided due to
        long range Coulomb interation.   See comments in the function.
        Here F is force within solvents particles.

        In RecomputeInitialFFTs(), pairwise long-range interactions in
        BHD.u2[][] are computed by ComputeFFTfromCoulomb().

        According to Page. 115 - 116:

                  0
          g(x) = g (x) exp[-u(x)]

        with

           0               LJ       Cs       Cl
          g (x) = exp[-β (V  (x) + V  (x) + V  (x)]

        then g(x) rewritten as:

                 ~ 0            Cl
          g(x) = g  (x) exp[-β V  (x) - u(x)]

        with

          ~ 0               LJ       Cs
          g  (x) = exp[-β (V  (x) + V  (x))]

        then BGY equation written as:

           ~              Cl
          Δu = K(g) + β ΔV

        with
                 ~ 0            Cl              ~ 0         ~
          g(x) = g  (x) exp[-β V  (x) - u(x)] = g  (x) exp[-u(x)]

                         ~
        the solution  of u can  be repsented by a  difference of
        two functions:

          ~   -    *
          u = u - u

        while:

           -              Cl       -
          Δu = K(g) + β ΔV   in Ω, u(∂Ω) = f,

            *           *
          Δu = 0 in Ω, u (∂Ω) = f,

        so after solving:

          Δu = K(g)

        we get:

          -          Cl
          u = u + β V

         Cl         Cl
        V      and V      are needed to calculated beforehand
         (A, M)     (B, M)

        and  sum to  the solution  by the  end of  each  iteration for
        solvent  site A  and B,  in  the code  hereafter.  These  site
        specific potentials can be  obtained by multiplying the common
        long-range Coulomb potential (stored in Vec uc) by the solvent
        site charge.
      */

      RecomputeInitialFFTs(&BHD, (damp > 0.0 ? damp : 0.0), 1.0);

      /* FIXME: avoid storing vectors  f[sl]_g2_fft on the grid across
         iterations. Only a scalar kernel is needed: */
      kernel (BHD.dc, BHD.PD, BHD.fs_g2_fft[0][0], NULL, ker_fft_S[0][0]);
      kernel (BHD.dc, BHD.PD, BHD.fs_g2_fft[0][1], NULL, ker_fft_S[0][1]); /* == [1][0] */
      kernel (BHD.dc, BHD.PD, BHD.fs_g2_fft[1][1], NULL, ker_fft_S[1][1]);

      kernel (BHD.dc, BHD.PD, BHD.fl_g2_fft[0][0], BHD.u2_fft[0][0], ker_fft_L[0][0]);
      kernel (BHD.dc, BHD.PD, BHD.fl_g2_fft[0][1], BHD.u2_fft[0][1], ker_fft_L[0][1]); /* == [1][0] */
      kernel (BHD.dc, BHD.PD, BHD.fl_g2_fft[1][1], BHD.u2_fft[1][1], ker_fft_L[1][1]);

      /* FIXME: what is  the point to split the  kernel in two pieces?
         Redefine S := S + L and forget about L */
      for (int i = 0; i < 2; i++)
        for (int j = 0; j <= i; j++) /* S := damp * L + damp_LJ * S */
          VecAXPBY (ker_fft_S[i][j], damp, damp_LJ, ker_fft_L[i][j]);

      /* Fill u0[0], u0[1] (alias BHD.g_ini[], also see the definition
         above) and  uc with VM_Coulomb_long.  No other  fields of the
         struct State except those passed explicitly are modified: */
      bgy3d_solute_field (&BHD,
                          2, solvent,
                          u0, uc, /* intent(out) */
                          n, solute,
                          density, /* void (*density)(...) */
                          (damp > 0.0 ? damp : 0.0), 1.0);

      for (int i = 0; i < 2; i++)
        {
          /*
            Historically  short-range  potential   is  scaled  by  the
            inverse  temperature and  the  code operates  with u(x)  =
            βv(x) insead of potential itself at many places:
          */
          VecScale (u0[i], beta);

          /*
            See pp.  116-177 in  thesis: boundary conditions (5.107) -
            (5.110):  first  impose   boundary  condition  then  solve
            laplacian equation  and substrate  from u0.  State  BHD is
            not   modified  by  these   calls.   Vec   work,  formally
            intent(out) in  these calls, is,  well, a work  array. Its
            value is ignored.
          */
          bgy3d_impose_laplace_boundary (&BHD, u0[i], work, BHD.x_lapl[i]);

          /* g :=  exp[-(u0 + du)] = g0 * exp(-du) */
          ComputeH2O_g (g[i], u0[i], du[i]);
        }

      real du_norm_old[2] = {0.0, 0.0}; /* Not sure if 0.0 as inital
                                           value is right.  */

      real a1 = a0;             /* loop-local variable */
      for (int iter = 0, mycount = 0, upwards = 0; iter < max_iter; iter++)
        {
          /* Every tenth iteration, starting with iter == 0: */
          const bool tenth = !(iter % 10);

          /*
            "a =  a1" is taken in  iteration 0, 10, 20,  etc.  "a1" is
            modified during the loop.

            "a =  a0" is  taken in iterations  1-9, 11-19,  etc.  "a0"
            remains unchanged during the loop.

            Note that in the first iteration a1 == a0.
          */
          const real a = tenth? a1 : a0;

          /*
            Some   functions,  such   as  bgy3d_solve_normalization(),
            Compute_dg_H2O_intra_ln(), and Compute_H2O_interS/_C() use
            preallocated fftw_complex arrays in State BHD for work but
            do not  re-define any  of the Vecs  in that  struct except
            those passed explicitly.

            Same     holds    for    Solve_NormalizationH2O_smallII(),
            bgy3d_impose_laplace_boundary()   to   the   best  of   my
            (limited) knowledge.
          */

          /* Compute FFT of g[] for all sites: */
          for (int i = 0; i < 2; i++)
            MatMult (BHD.fft_mat, g[i], g_fft[i]);

          /* for H, O in that order ... */
          for (int i = 0; i < 2; i++)
            {
              /*
                ... sum over H, O  in that order.  LJ, short- and long
                range  Coulomb, and  a so  called strange  addition is
                accounted for in the kernel.  First clear accumulator:
              */
              VecSet (du_acc_fft, 0.0);

              for (int j = 0; j < 2; j++) /* This increments the accumulator: */
                apply (BHD.dc, ker_fft_S[i][j], g_fft[j], beta * BHD.rhos[j], du_acc_fft);

              /*
                Compute IFFT of du_acc_fft for the current site. Other
                contributions  are  added  to  the real  space  du_acc
                below:
              */
              MatMultTranspose (BHD.fft_mat, du_acc_fft, du_acc);

              /*
                In the following the sum is  over all sites j /= i. It
                is supposed  to account for  intra-molecular site-site
                correlation  of the  solvent  sites due  to the  rigid
                bonds.  See e.g. Eq. (4.114), Jager Diss:
              */
              for (int j = 0; j < 2; j++)
                {
                  if (j == i) continue;
                  /*
                    This  body  is executed  exactly  once for  2-site
                    models.   FIXME:  in   a  more  general  case  the
                    distance  should be  the  intra-molecular distance
                    between sites i and j.

                    The  first   step  is  to   compute  normalization
                    functions ĝ(x), Eqs. (4.105), (4.106) and (4.110).
                    This already involves  a convolution integral with
                    the geometric factor δ(|x| - r) / 4πr².

                    Vec t_vec, is intent  (out) here and is re-used as
                    a work  vector to  pass the result  further.  Pass
                    the g(k), not g(x). Does one FFT^-1.
                  */
                  bgy3d_solve_normalization (&BHD, g_fft[i], r_HO, g[j], t_vec);

                  /*
                    The   next  step  is   to  use   the  "conditional
                    distribution" ĝ(x) again in a convolution integral
                    with the same geometric  factor δ(|x| - r) / 4πr²,
                    Eq. (4.114).

                    Vec t_vec is intent(in)  and work is intent (out).
                    Does one  FFT and one FFT^-1. Here  it is actually
                    possible to use the same Vec as in- and output:
                  */
                  Compute_dg_H2O_intra_ln (&BHD, t_vec, r_HO, work);

                  /* Add the contrinution of site j /= i to du[i]: */
                  VecAXPY (du_acc, 1.0, work);
                }

              /*
                Add Coulomb field uc scaled  by the site charge to the
                accumulator:
              */
              VecAXPY (du_acc, solvent[i].charge, uc);

              /*
                Vec   du_acc  and   BHD.x_lapl[i]   are  intent(inout)
                here. Try  to preserve  the values of  x_lapl[] across
                iterations to save time  in the iterative solver.  Vec
                t_vec is used as a work array.
              */
              bgy3d_impose_laplace_boundary (&BHD, du_acc, t_vec, BHD.x_lapl[i]);

              /*
                Mix du and du_new with a fixed ration "a":

                  du' = a * du_new + (1 - a) * du
                  norm = |du_new - du|
               */
              du_norm[i] = mix (du[i], du_acc, a, work); /* last arg is a temp */
            } /* over sites i */

          /*
            Now that du[] has been  computed using g[] of the previous
            iteration  one  can  safely  update  g[].   Compute  g  :=
            exp[-(u0 + du)], with a sanity check:
          */
          for (int i = 0; i < 2; i++)
            ComputeH2O_g (g[i], u0[i], du[i]);

          /*
           * Fancy step size control.
           *
           * FIXME: weired  logic.  The  code below appears  to modify
           * "upwards",  "a1", and  "mycount" by  eventually resetting
           * the latter to  zero. The goal might be  to tweak the real
           * coefficient  "a1"   depending  on  iteration   count  and
           * convergence.    Everytime  "mycount"   becomes   >20  the
           * coefficient "a1" is  changed. The compiler is complaining
           * that upwards  and du_nomr_old[] maybe  used uninitialized
           * here. At the moment I  am not able to confirm/reject that
           * claim.
           */

          /* Infinity (max-abs) norm of du[] over all site indices: */
          real norm8 = maxval (2, du_norm);

          mycount++;

          if ((iter - 1) % 10 && (du_norm_old[0] < du_norm[0] ||
                                  du_norm_old[1] < du_norm[1]))
            upwards = 1;
          else if (iter > 20 && !((iter - 1) % 10) && upwards == 0 &&
                  (du_norm_old[0] < du_norm[0] || du_norm_old[1] < du_norm[1]))
            {
              a1 /= 2.0;
              if (a1 < a0)
                a1 = a0;
              mycount = 0;
            }
          else
            upwards = 0;

          if (mycount > 20)
            {
              /* Scale the  coefficient "a1" up by a  factor, but make
                 sure it is not above 1.0. Reset mycount. */
              if (a1 <= 0.5)
                a1 *= 2.0;
              else
                a1 = 1.0;
              mycount = 0;
            }
          /* otherwise leave "a1" and "mycount" unchanged */

          for (int i = 0; i < 2; i++)
            du_norm_old[i] = du_norm[i];

          PetscPrintf (PETSC_COMM_WORLD, "%03d ", iter + 1);
          PetscPrintf (PETSC_COMM_WORLD, "a=%f ", a);
          PetscPrintf (PETSC_COMM_WORLD, "H=%e ", du_norm[0]);
          PetscPrintf (PETSC_COMM_WORLD, "O=%e ", du_norm[1]);

          /* Last argument to ComputeCharge() is a work array: */
          PetscPrintf (PETSC_COMM_WORLD, "Q=% e ",
                       ComputeCharge (PD, 2, solvent, g, work));

          PetscPrintf (PETSC_COMM_WORLD, "count=%3d upwards=%1d",
                       mycount, upwards);
          PetscPrintf (PETSC_COMM_WORLD, "\n");

          /* Exit  when any  of  du[]  does not  change  by more  than
             norm_tol: */
          if (norm8 <= norm_tol)
            {
              PetscPrintf (PETSC_COMM_WORLD,
                           "norm %e <= %e (norm-tol) in iteration %d < %d (max-iter)\n",
                           norm8, norm_tol, iter + 1, max_iter);
              break;
            }

        } /* iter loop */
      /*************************************/

      /* FIXME:  Debug  output  from  every iteration  with  different
         overall  scale  factors  damp/damp_LJ.  Remove when  no  more
         needed. */
      namecount++;
      sprintf(nameH, "vec0-%d.m", namecount-1);
      sprintf(nameO, "vec1-%d.m", namecount-1);

      PetscPrintf (PETSC_COMM_WORLD, "Writing files...");
      bgy3d_save_vec_ascii (nameH, g[0]); /* g_H */
      bgy3d_save_vec_ascii (nameO, g[1]); /* g_O */
      PetscPrintf (PETSC_COMM_WORLD, "done.\n");
      /************************************/

      /* Test for solvent electrostatic potential: */
      bgy3d_solvent_field (&BHD, 2, solvent, g, ve);

      /* Save du to binary file. FIXME: Why du and not g? */
      if (bgy3d_getopt_test ("--save-H2O"))
        {
          PetscPrintf (PETSC_COMM_WORLD, "Writing binary files...");
          bgy3d_save_vec ("du0.bin", du[0]); /* duH */
          bgy3d_save_vec ("du1.bin", du[1]); /* duO */
          PetscPrintf (PETSC_COMM_WORLD, "done.\n");
        }
    } /* damp loop */

  /* Clean up and exit ... */
  VecDestroy (du_acc_fft);

  for (int i = 0; i < 2; i++)
    {
      /* Delegated to the caller:
         VecDestroy (g[i]); */

      VecDestroy (du[i]);
      VecDestroy (g_fft[i]);

      for (int j = 0; j <= i; j++)
        {
          VecDestroy (ker_fft_S[i][j]);
          VecDestroy (ker_fft_L[i][j]);
        }
    }

  VecDestroy (du_acc);
  VecDestroy (work);
  VecDestroy (t_vec);
  VecDestroy (uc);
  VecDestroy (ve);

  finalize_state (&BHD);
}

/* This one emulates historical solver interface: */
Vec BGY3dM_solve_H2O_2site (const ProblemData *PD, Vec g_ini)
{
  (void) g_ini;                 /* FIXME: interface obligation */

  int n;                        /* number of solute sites */
  const Site *sites;            /* [n], array of sites */
  Vec g[2];                     /* solution */
  char name[200] = "hydrogen chloride"; /* default solute */

  /* Solutes name, HCl by default: */
  bgy3d_getopt_string ("--solute", name, sizeof(name));

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solute is %s.\n", name);

  /* Get the solute from the tables: */
  bgy3d_solute_get (name, &n, &sites);

  /* This does the  real work. Vec g[2] is  intent(out) in all senses,
     dont  forget   to  destroy   them.  Here  no   additional  charge
     distribution, so supply NULL for the function pointer: */
  bgy3d_solve_with_solute (PD, n, sites, NULL, g);

  /* Save final distribution, use binary format: */
  bgy3d_save_vec ("g0.bin", g[0]); /* gH */
  bgy3d_save_vec ("g1.bin", g[1]); /* gO */

  for (int i = 0; i < 2; i++)
    VecDestroy (g[i]);

  return PETSC_NULL;            /* fake, interface obligation */
}

Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini)
{
  (void) g_ini;                 /* FIXME: interface obligation */

  real a1, a, damp, damp_LJ;
  real count = 0.0;
  int iter;
  Vec g0H, g0O, dgH, dgO,  dg_new, dg_new2;
  Vec work;
  Vec g[2];
  Vec tH, tO, dg_newH, dg_newO;
  Vec uc;                       /* common for all sites */
  PetscScalar dgH_norm, dgO_norm;
  real dgH_old, dgO_old;
  int mycount=0, upwards, namecount=0;
  char nameH[20], nameO[20];

  Vec dg_histO, dg_histH;

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (3-site) equation ...\n");

  State BHD = initialize_state (PD);
  BHD.rhos[0] = 2.*BHD.rhos[0];

  if (r_HH < 0.0) {
    PetscPrintf(PETSC_COMM_WORLD,"Solvent not a 3-Site model!\n");
    exit(1);
  }

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Mixing parameter: */
  real a0 = PD->lambda;         /* not const */

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

#ifdef L_BOUNDARY
  /* Assemble Laplacian matrix and create KSP environment: */
  bgy3d_laplace_create (BHD.da, BHD.PD, &BHD.M, &BHD.ksp);
#endif
#ifdef L_BOUNDARY_MG
  /* Create KSP environment */
  InitializeDMMGSolver(&BHD);
#endif

  DACreateGlobalVector(BHD.da, &g[0]);
  DACreateGlobalVector(BHD.da, &g[1]);
  DACreateGlobalVector(BHD.da, &dgH);
  DACreateGlobalVector(BHD.da, &dgO);
  DACreateGlobalVector(BHD.da, &dg_new);
  DACreateGlobalVector(BHD.da, &dg_new2);
  DACreateGlobalVector(BHD.da, &work);

  DACreateGlobalVector(BHD.da, &tH);
  DACreateGlobalVector(BHD.da, &tO);

  DACreateGlobalVector(BHD.da, &dg_newH);
  DACreateGlobalVector(BHD.da, &dg_newO);

  DACreateGlobalVector(BHD.da, &uc); /* common for all sites */

  DACreateGlobalVector(BHD.da, &dg_histO);
  DACreateGlobalVector(BHD.da, &dg_histH);
  VecSet(dg_histH, 0.0);
  VecSet(dg_histO, 0.0);

  g0H=BHD.g_ini[0];
  g0O=BHD.g_ini[1];

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dg_new,0.0);

  if (bgy3d_getopt_test ("--from-g2")) {
      ComputedgFromg (dgH, g0H, BHD.g2[0][1]);
      ComputedgFromg (dgO, g0O, BHD.g2[1][1]);
  }

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-H2O")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      dgH = bgy3d_load_vec ("dg0.bin"); /* dgH */
      dgO = bgy3d_load_vec ("dg1.bin"); /* dgO */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  }

  for( damp=damp_start; damp <=1; damp+=0.1)
    {
      if (damp <= -0.01)
        {
          damp_LJ = 0.0;
        }
      else if (damp==0.0)
        {
          damp_LJ = 1.0;
        }
      else
        {
          damp_LJ = 1.0;
          count += 1.0;
          a0 = 0.1 / (count + 5.0);
        }

      RecomputeInitialFFTs(&BHD, (damp > 0.0 ? damp : 0.0), damp_LJ);

      /* Compute solute field in this block:*/
      {
        int n;                        /* number of solute sites */
        const Site *sites;            /* [n], array of sites */

        /* Get the solute from the tables: */
        bgy3d_solute_get ("butanoic acid", &n, &sites);

        /* This does the real work: */
        bgy3d_solute_field (&BHD,
                            2, solvent,
                            BHD.g_ini, uc, /* intent(out) */
                            n, sites,
                            NULL, /* void (*density)(...) */
                            (damp > 0.0 ? damp : 0.0), damp_LJ);
      }

      PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);

      /* Historically   short-range  potential   is   stored  with   a
         factor: */
      for (int i = 0; i < 2; i++)
        VecScale (BHD.g_ini[i], BHD.PD->beta);

      bgy3d_impose_laplace_boundary (&BHD, g0H, tH, BHD.x_lapl[0]);
      bgy3d_impose_laplace_boundary (&BHD, g0O, tH, BHD.x_lapl[1]);

      /* g=g0*exp(-dg) */
      ComputeH2O_g( g[0], g0H, dgH);
      ComputeH2O_g( g[1], g0O, dgO);

      a=a0;
      a1=a0;
      for(iter=0; iter<max_iter; iter++)
        {

          if( !(iter%10) && iter>0 )
            a=a1;
          else
            a=a0;

          PetscPrintf (PETSC_COMM_WORLD, "%03d: ", iter + 1);


          /* H */
          VecSet(dg_new,0.0);
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[0][1], g[1], BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[0][0], g[0], BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][1], g[1], BHD.u2_fft[0][1], damp * BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][0], g[0], BHD.u2_fft[0][0], damp * BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_smallII (&BHD, g[0], r_HH, g[0], tH , dg_new2, work);

          Compute_dg_H2O_intra_ln(&BHD, tH, r_HH, dg_new2);
          VecCopy (dg_new2, work); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);
          Solve_NormalizationH2O_smallII (&BHD, g[0], r_HO, g[1], tO , dg_new2, work);

          Compute_dg_H2O_intra_ln(&BHD, tO, r_HO, dg_new2);
          VecCopy (dg_new2, work); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

          /* Copy  the electrostatic potential  scaled by  the solvent
             site charges into predefined locations: */
          VecAXPY(dg_new, solvent[0].charge, uc);

          bgy3d_impose_laplace_boundary (&BHD, dg_new, tH, BHD.x_lapl[0]);

          VecCopy(dg_new, dg_newH);

          /* O */
          VecSet(dg_new,0.0);
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[1][1], g[1], BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[0][1], g[0], BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[1][1], g[1], BHD.u2_fft[1][1], damp * BHD.rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][1], g[0], BHD.u2_fft[0][1], damp * BHD.rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_smallII (&BHD, g[1], r_HO, g[0], tH , dg_new2, work);
          Compute_dg_H2O_intra_ln(&BHD, tH, r_HO, dg_new2);
          VecCopy (dg_new2, work); /* FIXME: need that? */

          VecAXPY(dg_new, 2.0, dg_new2);

          /* Copy  the electrostatic potential  scaled by  the solvent
             site charges into predefined locations: */
          VecAXPY(dg_new, solvent[1].charge, uc);

          bgy3d_impose_laplace_boundary (&BHD, dg_new, tH, BHD.x_lapl[1]);

          VecCopy(dg_new, dg_newO);

          /* Move dgH */
          VecCopy(dgH, work);
          VecAXPBY(dgH, a, (1-a), dg_newH);
          VecAXPY(work, -1.0, dgH);
          VecNorm(work, NORM_INFINITY, &dgH_norm);

          PetscPrintf(PETSC_COMM_WORLD,"H= %e (a=%f) ", dgH_norm/a, a);

          /* Move dgO */
          if (1)
            {
              VecCopy(dgO, work);
              VecAXPBY(dgO, a, (1-a), dg_newO);
              VecAXPY(work, -1.0,  dgO);
              VecNorm(work, NORM_INFINITY, &dgO_norm);
              PetscPrintf(PETSC_COMM_WORLD,"O= %e (a=%f) ", dgO_norm/a, a);
            }
          ComputeH2O_g( g[0], g0H, dgH);
          ComputeH2O_g( g[1], g0O, dgO);

          /* Last argument to ComputeCharge() is a work array: */
          PetscPrintf (PETSC_COMM_WORLD, "Q=% e ",
                       ComputeCharge (PD, 2, solvent, g, work));

          /* (fancy) step size control */
          mycount++;
          if( ((iter-1)%10) &&
              (dgH_old<dgH_norm/a || dgO_old<dgO_norm/a ) )
            {
              upwards=1;
            }
          else if(iter>20 && !((iter-1)%10) && upwards==0 &&
                  (dgH_old<dgH_norm/a || dgO_old<dgO_norm/a ) )
            {
              a1 /=2.;
              if(a1<a0)
                a1=a0;
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
          dgH_old = dgH_norm/a;
          dgO_old = dgO_norm/a;

          PetscPrintf(PETSC_COMM_WORLD,"\n");

          if(dgH_norm/a<=norm_tol &&  dgO_norm/a<=norm_tol)
            break;
        }

      /* output */
      namecount++;
      sprintf(nameH, "vec0-%d.m", namecount-1);
      sprintf(nameO, "vec1-%d.m", namecount-1);

      PetscPrintf(PETSC_COMM_WORLD,"Writing files...");
      bgy3d_save_vec_ascii (nameH, g[0]); /* g_H */
      bgy3d_save_vec_ascii (nameO, g[1]); /* g_O */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");

      /* save g to binary file */
      if (bgy3d_getopt_test ("--save-H2O")) {
          PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
          bgy3d_save_vec ("dg0.bin", dgH); /* dgH */
          bgy3d_save_vec ("dg1.bin", dgH); /* dgO */
          PetscPrintf(PETSC_COMM_WORLD,"done.\n");
      }
    }

  VecDestroy(g[0]);
  VecDestroy(g[1]);
  VecDestroy(dgH);
  VecDestroy(dgO);
  VecDestroy(dg_new);
  VecDestroy(dg_new2);
  VecDestroy(work);

  VecDestroy(tH);
  VecDestroy(tO);

  VecDestroy (uc);

  VecDestroy(dg_newH);
  VecDestroy(dg_newO);
  VecDestroy(dg_histH);
  VecDestroy(dg_histO);

  finalize_state (&BHD);

  return PETSC_NULL;
}
