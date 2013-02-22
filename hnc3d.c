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

typedef struct HNC3dData
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
} HNC3dData;


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


HNC3dData* HNC3dData_malloc(const ProblemData *PD)
{
  HNC3dData *HD = malloc (sizeof (*HD));

  HD->PD = PD;

  /* Initialize  parallel  stuff,  fftw  +  petsc.  Data  distribution
     depends on the grid dimensions N[] and number of processors.  All
     other arguments are intent(out): */
  bgy3d_fft_mat_create (PD->N, &HD->fft_mat, &HD->da, &HD->dc);

  DA da = HD->da;

  /* Create global vectors */
  DACreateGlobalVector (da, &HD->pot);

  /* FIXME:   this  is   abused  to   get  both   solvent-solvent  and
     solute-solvent interactions: */
  solute_field (HD->da, HD->PD, HD->pot);

  return HD;
}


static void HNC3dData_free (HNC3dData *HD)
{
  VecDestroy (HD->pot);

  DADestroy (HD->da);
  DADestroy (HD->dc);
  MatDestroy (HD->fft_mat);

  free(HD);
}


/*
  A function  that takes a  context, an input  Vec x and  computes the
  residual Vec r.  An example is the "error  vector" of the non-linear
  equation is  the difference  between the input  and the output  of a
  "fixpoint" iteration as a function of input:

    r = x    - x
         out    in
*/
typedef void (*Function) (void *ctx, Vec x, Vec r);


/*
  For solving HNC equation with Newton. Except of x everything else in
  the  closure context is  considered read  only input  or intemediate
  terms depending  on x.  That  is when looking for  total correlation
  function h,  the direct correlation  function should be fixed  (or a
  function of h).
*/
static void snes_solve (void *ctx, Function F, Vec x)
{
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

  /* SNES needs a place to store residual: */
  Vec r;
  VecDuplicate (x, &r);

  /* SNES functions should obey this interface: */
  PetscErrorCode F1 (SNES snes, Vec x, Vec r, void *ctx)
  {
    (void) snes;                /* unused */
    F (ctx, x, r);              /* assumes ctx is a Context* */
    return 0;
  }
  SNESSetFunction (snes, r, F1, ctx); /* Pass Context* as ctx */

  /*
    Runtime  options will  override default  parameters.   FIXME: note
    that the  call to SNESSetJacobian()  is missing here.   It appears
    that  one has to  request a  "matrix-free" approximation  from the
    command line  with "-snes_mf". Otherwise the  next call terminates
    with an error message saying "Matrix must be set first"!
  */
  SNESSetFromOptions (snes);

  /* Solve  problem F(x)  = 0.  PETSC_NULL indicates  that the  rhs is
     0: */
  SNESSolve (snes, PETSC_NULL, x);

  /* Write  out  solution.   FIXME:   and  what  was  the  purpose  of
     SNESSolve()? */
  // SNESGetSolution (snes, &x);

  VecDestroy (r);

  SNESDestroy (snes);
}


static void picard_solve (const ProblemData *PD, void *ctx, Function F, Vec x)
{
  /* Mixing parameter */
  const real lambda = PD->lambda;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* Convergence threshold: */
  const real norm_tol = PD->norm_tol;

  /* A place to store residual: */
  Vec dx;
  VecDuplicate (x, &dx);

  /* Find an x such that dx as returned by F (ctx, x, dx) is zero: */
  for (int k = 0; k < max_iter; k++)
    {
      F (ctx, x, dx);

      /* Simple mixing: x = lambda * x + (1 - lambda) * x_old */
      VecAXPY (x, lambda, dx);

      const real norm = bgy3d_vec_norm (dx);

      PetscPrintf (PETSC_COMM_WORLD, "%03d: norm of difference: %e\t%f\n",
                   k + 1, norm, lambda);

      if (norm < norm_tol)
        break;
    }
  VecDestroy (dx);
}


/* h := exp (-β v) - 1 */
static void compute_h (real beta, Vec v, Vec h)
{
  real pure f (real v)
  {
    /* exp (-beta * v) - 1 */
    return expm1 (-beta * v);
  }
  bgy3d_vec_map1 (h, f, v);
}


/*
  Hypernetted   Chain  (HNC)  closure   relation  to   compute  direct
  correlation function c in real space.  See OZ equation below for the
  second  relation between  two unknowns  c and  γ. Other  sources use
  latin "t"  to denote that difference  and we will use  that to avoid
  confusion with distributions functions:

    c := exp (-β v + γ) - 1 - γ
*/
static void compute_c (real beta, Vec v, Vec t, Vec c)
{
  real pure f (real v, real t)
  {
    /* exp (-beta * v + t) - 1.0 - t */
    return expm1 (-beta * v + t) - t;
  }
  bgy3d_vec_map2 (c, f, v, t);
}


/*
  Use the k-representation of Ornstein-Zernike (OZ) equation

    h = c + ρ c * h

  to compute γ =  h - c form c:

          2
    γ = ρc  / (1 - ρc)

  If you scale  c by h3 beforehand or  pass rho' = rho *  h3 and scale
  the  result by h3  in addition,  you will  compute exactly  what the
  older version of the function did:
*/
static void compute_t (real rho, Vec c_fft, Vec t_fft)
{
  complex pure f (complex c)
  {
    return rho * (c * c) / (1.0 - rho * c);
  }
  bgy3d_vec_fft_map1 (t_fft, f, c_fft);
}


/*
  HNC iteration for a direct correlation where the total correlation
  is considered a function of that:

  c    ->  dc = c    - c
    in           out    in
*/
typedef struct Ctx_c
{
  HNC3dData *HD;
  Vec t;                        /* real */
  Vec c_fft, t_fft;             /* complex */
} Ctx_c;


static void iterate_c (Ctx_c *ctx, Vec c, Vec dc)
{
  const ProblemData *PD = ctx->HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real L = PD->interval[1] - PD->interval[0];
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  MatMult (ctx->HD->fft_mat, c, ctx->c_fft);

  VecScale (ctx->c_fft, h3);

  compute_t (rho, ctx->c_fft, ctx->t_fft);

  MatMultTranspose (ctx->HD->fft_mat, ctx->t_fft, ctx->t);

  VecScale (ctx->t, 1.0/L/L/L);

  compute_c (beta, ctx->HD->pot, ctx->t, dc);

  VecAXPY (dc, -1.0, c);
}


/* Solving  for  c  and  h(c)  of HNC  equation  with  Newton.  Direct
   correlation c appears as a primary variable here: */
Vec hnc3d_solve_newton (const ProblemData *PD, Vec g_ini)
{
  assert(g_ini==PETSC_NULL);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving 3d-HNC equation. Newton iteration.\n");

  HNC3dData *HD = HNC3dData_malloc (PD);
  Vec t = bgy3d_vec_create (HD->da);
  Vec t_fft = bgy3d_vec_create (HD->dc); /* complex */
  Vec c_fft = bgy3d_vec_create (HD->dc); /* complex */

  /* Create intial guess: */
  Vec c = bgy3d_vec_create (HD->da);
  VecSet (c, 0.0);

  /*
    Find a  c such that dc as  returned by iterate_c (&ctx,  c, dc) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: struct Ctx_c* vs. void*:
  */
  {
    struct Ctx_c ctx = {HD, t, t_fft, c_fft};
    snes_solve (&ctx, (Function) iterate_c, c);
  }

  /* g = γ + c + 1, store in Vec t: */
  VecAXPY (t, 1.0, c);
  VecShift (t, 1.0);

  bgy3d_vec_save ("c00.bin", c);
  bgy3d_vec_save ("g00.bin", t);

  /* free stuff */
  /* Delegated to the caller: VecDestroy (t); */
  VecDestroy (t_fft);
  VecDestroy (c);
  VecDestroy (c_fft);

  HNC3dData_free (HD);

  return t;
}


/* Solve h and c of HNC equation simultaneously, fixpoint iteration */
Vec hnc3d_solve_picard (const ProblemData *PD, Vec g_ini)
{
  assert(g_ini==PETSC_NULL);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving 3d-HNC equation. Fixpoint iteration.\n");

  HNC3dData *HD = HNC3dData_malloc (PD);
  Vec t = bgy3d_vec_create (HD->da);
  Vec t_fft = bgy3d_vec_create (HD->dc); /* complex */
  Vec c_fft = bgy3d_vec_create (HD->dc); /* complex */

  /* Create intial guess: */
  Vec c = bgy3d_vec_create (HD->da);
  Vec dc = bgy3d_vec_create (HD->da);
  VecSet (c, 0.0);

  /*
    Find a  c such that dc as  returned by iterate_c (&ctx,  c, dc) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: struct Ctx_c* vs. void*:
  */
  {
    /* Work area for iterate_c(): */
    struct Ctx_c ctx = {HD, t, t_fft, c_fft};
    picard_solve (PD, &ctx, (Function) iterate_c, c);
  }

  /* g = γ + c + 1, store in Vec t: */
  VecAXPY (t, 1.0, c);
  VecShift (t, 1.0);

  bgy3d_vec_save ("c00.bin", c);
  bgy3d_vec_save ("g00.bin", t);

  /* free stuff */
  /* Delegated to the caller: VecDestroy (t); */
  VecDestroy (t_fft);
  VecDestroy (c);
  VecDestroy (c_fft);
  VecDestroy (dc);

  HNC3dData_free (HD);

  return t;
}


/*
  HNC iteration for a fixed direct correlation:

  h    ->  dh = h    - h
    in           out    in
*/
typedef struct Ctx_h
{
  HNC3dData *HD;
  Vec t;                      /* real */
  Vec c_fft, h_fft, ch_fft;     /* complex */
} Ctx_h;

static void iterate_h (Ctx_h *ctx, Vec h, Vec dh)
{
  const ProblemData *PD = ctx->HD->PD;

  Vec t = ctx->t;           /* temp */
  Vec c_fft = ctx->c_fft;       /* fixed solvent kernel */
  Vec h_fft = ctx->h_fft;       /* temp */
  Vec ch_fft = ctx->ch_fft;     /* temp */

  const real rho = PD->rho;
  const real beta = PD->beta;

  /* fft(h) */
  MatMult (ctx->HD->fft_mat, h, h_fft);

  /* fft(h)*fft(c) */
  complex pure mul (complex x, complex y)
  {
    return x * y;
  }
  bgy3d_vec_fft_map2 (ch_fft, mul, c_fft, h_fft);

  /* v = fft^-1(fft(c)*fft(h)) */
  MatMultTranspose (ctx->HD->fft_mat, ch_fft, t);

  VecScale (t, PD->h[0] * PD->h[1] * PD->h[2] / PD->N[0] / PD->N[1] / PD->N[2]);

  /*
    The new candidate for the total correlation

    h = exp (-βU + ρv) - 1

    stored in Vec dh:
  */
  real pure h_out (real u, real g)
  {
    return expm1 (-beta * u + rho * g);
  }
  bgy3d_vec_map2 (dh, h_out, ctx->HD->pot, t);

  /*
    dh := h    - h
           out    in
  */
  VecAXPY (dh, -1.0, h);
}

static void solvent_kernel (HNC3dData *HD, Vec c_fft)
{
  Vec c = bgy3d_vec_create (HD->da);

  /* Load c_1d from file: */
  if (bgy3d_getopt_test ("--from-radial-g2")) /* FIXME: better name? */
    bgy3d_vec_read_radial (HD->da, HD->PD, "c1dfile", c);
  else
    bgy3d_vec_read ("c00.bin", c);

  MatMult (HD->fft_mat, c, c_fft);

  VecDestroy (c);
}


/* solving for h only of HNC equation with Newton */
/* c appears as an input here */
Vec hnc3d_solute_solve_newton (const ProblemData *PD, Vec g_ini)
{
  assert(g_ini==PETSC_NULL);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving 3d-HNC equation. Newton iteration. Fixed c.\n");

  HNC3dData *HD = HNC3dData_malloc (PD);
  Vec t = bgy3d_vec_create (HD->da);
  Vec c_fft = bgy3d_vec_create (HD->dc);  /* complex */
  Vec h_fft = bgy3d_vec_create (HD->dc);  /* complex */
  Vec ch_fft = bgy3d_vec_create (HD->dc); /* complex */

  /* Get the solvent-solvent direct correlation function: */
  solvent_kernel (HD, c_fft);

  /* Create global vectors */
  Vec h;
  VecDuplicate (HD->pot, &h);

  /* Set initial guess */
  compute_h (HD->PD->beta, HD->pot, h);

  /*
    Find  an h such  that dh  as returned  by iterate  (HD, h,  dh) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: Ctx_h* vs. void*:
  */
  {
    /* Work area for iterate_h(): */
    struct Ctx_h ctx =
      {
        .HD = HD,
        .t = t,
        .c_fft = c_fft,
        .h_fft = h_fft,
        .ch_fft = ch_fft
      };

    snes_solve (&ctx, (Function) iterate_h, h);
  }

  /* free stuff */
  /* Delegated to the caller: VecDestroy (h) */
  VecDestroy (t);
  VecDestroy (h_fft);
  VecDestroy (c_fft);
  VecDestroy (ch_fft);

  HNC3dData_free (HD);

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
               "Solving 3d-HNC equation. Fixpoint iteration. Fixed c.\n");

  HNC3dData *HD = HNC3dData_malloc (PD);
  Vec t = bgy3d_vec_create (HD->da);
  Vec c_fft = bgy3d_vec_create (HD->dc);  /* complex */
  Vec h_fft = bgy3d_vec_create (HD->dc);  /* complex */
  Vec ch_fft = bgy3d_vec_create (HD->dc); /* complex */

  /* Get the solvent-solvent direct correlation function: */
  solvent_kernel (HD, c_fft);

  Vec h, dh;
  VecDuplicate (HD->pot, &h);
  VecDuplicate (HD->pot, &dh);

  /* Set initial guess */
  compute_h (HD->PD->beta, HD->pot, h);

  /*
    Find a  c such that dc as  returned by iterate_c (&ctx,  c, dc) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: struct Ctx_c* vs. void*:
  */
  {
    /* Work area for iterate_h(): */
    struct Ctx_h ctx =
      {
        .HD = HD,
        .t = t,
        .c_fft = c_fft,
        .h_fft = h_fft,
        .ch_fft = ch_fft
      };

    picard_solve (PD, &ctx, (Function) iterate_h, h);
  }

  /* g = h + 1 */
  VecShift (h, 1.0);

  bgy3d_vec_save ("g0.bin", h);

  /* free stuff */
  /* Delegated to the caller: VecDestroy (h) */
  VecDestroy (t);
  VecDestroy (c_fft);
  VecDestroy (h_fft);
  VecDestroy (ch_fft);
  VecDestroy (dh);

  HNC3dData_free (HD);

  return h;
}
