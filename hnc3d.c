/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-vec.h"          /* bgy3d_vec_create() */
#include "bgy3d-force.h"        /* Lennard_Jones() */
#include "hnc3d.h"
#include <math.h>               /* expm1() */

typedef struct HNC3dDataStruct
{
  /*
    These are array  descriptors for real and complex  vectors and the
    FFT matrix.  The  data distribution of Petsc vectors  and FFTW MPI
    needs to be consistent, so  that these three should be constructed
    accordingly:
  */
  DA da, dc;
  Mat fft_mat;

  /* Immutable command line parameters are stored here: */
  const ProblemData *PD;

  Vec pot;
  Vec h_ini;

  /* things for arbitrary molecule shape */
  Vec c, v;
  Vec c_fft, h_fft, ch_fft;     /* complex */
} *HNC3dData;


static void solute_field (const DA da, const ProblemData *PD, Vec pot)
{
  real **x_M;
  int N_M;

  real interval[2];
  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  real h[3];
  FOR_DIM
    h[dim] = PD->h[dim];

  const real epsilon = 1.0;
  const real sigma = 1.0;

  /* Load molecule from file */
  x_M = Load_Molecule (&N_M);

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
          /* set force vector */
          /* loop over particles and grid */
          for (int k = 0; k < N_M; k++)
            {
              real r[3];
              FOR_DIM
                r[dim] = i[dim] * h[dim] + interval[0] - x_M[k][dim];

              const real r_s = sqrt (SQR (r[0]) + SQR (r[1]) + SQR (r[2]));

              pot_[i[2]][i[1]][i[0]] +=
                Lennard_Jones (r_s, epsilon, sigma);
            }
        }
  DAVecRestoreArray (da, pot, &pot_);

  Molecule_free (x_M, N_M);
}


HNC3dData HNC3dData_malloc(const ProblemData *PD)
{
  HNC3dData HD;
  DA da;

  HD = (HNC3dData) malloc(sizeof(*HD));

  HD->PD = PD;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &HD->fft_mat, &HD->da, &HD->dc);

  da = HD->da;

  /* Create global vectors */
  DACreateGlobalVector (da, &HD->pot);
  VecDuplicate (HD->pot, &HD->c);
  VecDuplicate (HD->pot, &HD->v);
  VecDuplicate (HD->pot, &HD->h_ini);

  /* FIXME:   this  is   abused  to   get  both   solvent-solvent  and
     solute-solvent interactions: */
  solute_field (HD->da, HD->PD, HD->pot);

  const real beta = PD->beta;

  real h0 (real v)
  {
    /* exp (-beta * v) - 1.0; */
    return expm1 (-beta * v);
  }
  bgy3d_vec_map1 (HD->h_ini, h0, HD->pot);

  HD->c_fft = bgy3d_vec_create (HD->dc);
  HD->h_fft = bgy3d_vec_create (HD->dc);
  HD->ch_fft = bgy3d_vec_create (HD->dc);

  return HD;
}


static void HNC3dData_free(HNC3dData HD)
{
  VecDestroy(HD->pot);
  VecDestroy(HD->c);
  VecDestroy(HD->v);
  DADestroy (HD->da);
  DADestroy (HD->dc);
  MatDestroy (HD->fft_mat);

  VecDestroy (HD->c_fft);
  VecDestroy (HD->h_fft);
  VecDestroy (HD->ch_fft);


  free(HD);
}


/* c := exp (-β v + γ) - 1 - γ */
static void Compute_c_HNC (real beta, Vec v, Vec g, Vec c)
{
  real pure f (real v, real g)
  {
    /* exp (-beta * v + g) - 1.0 - g */
    return expm1 (-beta * v + g) - g;
  }
  bgy3d_vec_map2 (c, f, v, g);
}


/*
  In k-representation compute

        2
  γ = ρc  / (1 - ρc)

  If you scale  c by h3 beforehand or  pass rho' = rho *  h3 and scale
  the  result by h3  in addition,  you will  compute exactly  what the
  older version of the function did:
*/
static void Compute_cgfft (real rho, Vec c_fft, Vec g_fft)
{
  complex pure f (complex c)
  {
    return rho * (c * c) / (1.0 - rho * c);
  }
  bgy3d_vec_fft_map1 (g_fft, f, c_fft);
}


/* Solve h and c of HNC equation simultaneously, fixpoint iteration */
Vec hnc3d_solve (const ProblemData *PD, Vec g_ini)
{
  HNC3dData HD;
  Vec c, c_old, g, g_old, gg;
  real g_norm, iL3;
  int k, n[3], x[3];
  Vec c_fft, cg_fft;

  assert(g_ini==PETSC_NULL);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving 3d-HNC equation. Fixpoint iteration.\n");

  /* Mixing parameter */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* Convergence threshold: */
  const real norm_tol = PD->norm_tol;

  HD = HNC3dData_malloc(PD);

  DAGetCorners (HD->da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  iL3 = 1./pow(PD->interval[1]-PD->interval[0],3);

  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  VecDuplicate(HD->pot, &c);
  VecDuplicate(HD->pot, &g);
  VecDuplicate(HD->pot, &g_old);
  VecDuplicate(HD->pot, &c_old);
  VecDuplicate(HD->pot, &gg);

  /* Set initial guess */
  VecSet(g,0.0);
  VecSet(c,0.0);

  /* set fft data */
  c_fft = bgy3d_vec_create (HD->dc);
  cg_fft = bgy3d_vec_create (HD->dc);

  /* do the iteration */
  for(k=0; k<max_iter; k++)
    {
      /* if(k>3) */
      /*   lambda =0.1; */

      /* set g_old=g */
      VecCopy(g, g_old);
      VecCopy(c, c_old);

      Compute_c_HNC(HD->PD->beta, HD->pot, g, c);

      /* simple mixing: c = lambda*c+(1-lambda)*c_old */
      VecAXPBY(c, (1-lambda), lambda, c_old);

      MatMult (HD->fft_mat, c, c_fft);

      VecScale (c_fft, h3);

      Compute_cgfft(PD->rho, c_fft, cg_fft);

      MatMultTranspose (HD->fft_mat, cg_fft, g);

      VecScale(g, iL3);
/*       VecView(c,PETSC_VIEWER_STDERR_WORLD);   */
/*       exit(1);   */
      /* gg=g-g_old */
      VecWAXPY(gg, -1.0, g_old, g);

      VecNorm(gg, NORM_2, &g_norm);
      PetscPrintf (PETSC_COMM_WORLD, "%03d: norm of difference: %e\t%f\n",
                   k + 1, g_norm, lambda);

      if (g_norm < norm_tol)
        break;
    }

  /* g= gamma+c+1 */
  VecAXPY(g, 1.0, c);
  VecShift(g, 1.0);

  //VecCopy(c,g);

  bgy3d_vec_save ("c00.bin", c);
  bgy3d_vec_save ("g00.bin", g);

  /* free stuff */
  VecDestroy(c);
  VecDestroy(gg);
  VecDestroy(c_old);
  VecDestroy(g_old);

  VecDestroy (c_fft);
  VecDestroy (cg_fft);

  HNC3dData_free(HD);

  return g;
}

/*
  HNC iteration for a fixed direct correlation:

  h    ->  dh = h    - h
    in           out    in
 */
static void iterate (HNC3dData HD, Vec h, Vec dh)
{
  const ProblemData *PD = HD->PD;

  /* fft(h) */
  MatMult (HD->fft_mat, h, HD->h_fft);

  Vec c_fft = HD->c_fft;
  Vec h_fft = HD->h_fft;
  Vec ch_fft= HD->ch_fft;

  const real rho = HD->PD->rho;
  const real beta= HD->PD->beta;

  /* fft(h)*fft(c) */
  complex pure mul (complex x, complex y)
  {
    return x * y;
  }
  bgy3d_vec_fft_map2 (ch_fft, mul, c_fft, h_fft);

  /* v = fft^-1(fft(c)*fft(h)) */
  MatMultTranspose (HD->fft_mat, ch_fft, HD->v);

  VecScale (HD->v, PD->h[0] * PD->h[1] * PD->h[2] / PD->N[0] / PD->N[1] / PD->N[2]);

  /*
    The new candidate for the total correlation

    h = exp (-βU + ρv) - 1

    stored in Vec dh:
  */
  real pure h_out (real u, real v)
  {
    return expm1 (-beta * u + rho * v);
  }
  bgy3d_vec_map2 (dh, h_out, HD->pot, HD->v);

  /*
    dh := h    - h
           out    in
  */
  VecAXPY (dh, -1.0, h);
}

/* F for solving HNC equation with Newton */
/* c appears as an input here */
static PetscErrorCode snes_function (SNES snes, Vec h, Vec f, void *pa)
{
  (void) snes;                  /* interface obligation */

  HNC3dData HD = (HNC3dData) pa;

  /*
    The "error vector" of the non-linear equation is:

    f := h    - h
          out    in
  */
  iterate (HD, h, f);

  return 0;                     /* interface obligation */
}


static void solvent_kernel (HNC3dData HD, Vec c, Vec c_fft)
{
  /* Load c_1d from file: */
  if (bgy3d_getopt_test ("--from-radial-g2")) /* FIXME: better name? */
    bgy3d_vec_read_radial (HD->da, HD->PD, "c1dfile", c);
  else
    bgy3d_vec_read ("c00.bin", c);

  MatMult (HD->fft_mat, c, c_fft);
}


/* solving for h only of HNC equation with Newton */
/* c appears as an input here */
Vec hnc3d_solute_solve_newton(const ProblemData *PD, Vec g_ini)
{
  assert(g_ini==PETSC_NULL);

  HNC3dData HD = HNC3dData_malloc(PD);

  /* Get the solvent-solvent direct correlation function: */
  solvent_kernel (HD, HD->c, HD->c_fft);

  /* Create global vectors */
  Vec F, h;
  VecDuplicate (HD->pot, &F);
  VecDuplicate (HD->pot, &h);

  /* Initial  guess. FIXME: this  is only  meaningfull for  solvent as
     solute: */
  VecCopy (HD->c, h);
  // VecSet (h, 0.0);

  /* Create the snes environment */
  SNES snes;
  SNESCreate (PETSC_COMM_WORLD, &snes);

  KSP ksp;
  SNESGetKSP (snes, &ksp);

  PC pc;
  KSPGetPC (ksp, &pc);

  /* set rtol, atol, dtol, maxits */
  KSPSetTolerances (ksp, 1.0e-5, 1.0e-50, 1.0e+5, 1000);

  /* line search: SNESLS, trust region: SNESTR */
  SNESSetType (snes, SNESLS);

  /* set preconditioner: PCLU, PCNONE, PCJACOBI... */
  PCSetType (pc, PCNONE);

/*   ComputeHNC_F(snes, g, F, (void*) HD);  */
/*   VecView(F,PETSC_VIEWER_STDERR_WORLD);  */
/*   exit(1);   */

  SNESSetFunction (snes, F, snes_function, HD);

  /*
    Runtime  options will  override default  parameters.   FIXME: note
    that the  call to SNESSetJacobian()  is missing here.   It appears
    tha  one has  to request  a "matrix-free"  approximation  from the
    command line  with "-snes_mf". Otherwise the  next call terminates
    with an error message saying "Matrix must be set first"!
  */
  SNESSetFromOptions (snes);

  /* solve problem */
  SNESSolve (snes, PETSC_NULL, h);

  /* write out solution */
  SNESGetSolution (snes, &h);

  /* free stuff */
  HNC3dData_free (HD);

  VecDestroy (F);

  SNESDestroy (snes);

  /* g := h + 1 */
  VecShift (h, 1.0);

  bgy3d_vec_save ("g0.bin", h);

  return h;

}


/* Solving for h of HNC eqauation with fixpoint iteration */
/* c is input */
Vec hnc3d_solute_solve_picard (const ProblemData *PD, Vec g_ini)
{
  assert (g_ini == PETSC_NULL);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving 3d-HNC equation by fixpoint iteration. Fixed c.\n");

  /* Mixing parameter: */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  HNC3dData HD = HNC3dData_malloc(PD);

  /* Get the solvent-solvent direct correlation function: */
  solvent_kernel (HD, HD->c, HD->c_fft);

  Vec h, dh;
  VecDuplicate (HD->pot, &h);
  VecDuplicate (HD->pot, &dh);

  /* Set initial guess */
  VecCopy (HD->h_ini, h);
  //VecSet(h,0.0);

  /* do the iteration */
  for (int k = 0; k < max_iter; k++)
    {
      /* dh := h_out - h_in */
      iterate (HD, h, dh);

      /* Simple mixing: h = lambda * h_out + (1 - lambda) * h_in */
      VecAXPY (h, lambda, dh);

      real norm;
      VecNorm (dh, NORM_INFINITY, &norm);
      PetscPrintf (PETSC_COMM_WORLD, "%03d: norm of difference: %e\t%f\n",
		   k + 1, norm, lambda);

      if (norm < norm_tol)
        break;
    }

  /* g = h + 1 */
  VecShift (h, 1.0);

  bgy3d_vec_save ("g0.bin", h);

  /* free stuff */
  /* Delegated to the caller: VecDestroy (h); */
  VecDestroy (dh);

  HNC3dData_free (HD);

  return h;
}
