/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-vec.h"          /* bgy3d_vec_create() */
#include "bgy3d-solutes.h"      /* struct Site */
#define G 1.2
#include "bgy3d-force.h"        /* Lennard_Jones() */
#include "bgy3d-snes.h"         /* bgy3d_snes_newton() */
#include "hnc3d-sles.h"         /* hnc3d_sles_zgesv() */
#include "hnc3d.h"
#include <math.h>               /* expm1() */

// #define HNC3D_T                 /* use γ as primary variable */


static void pair (const DA da, const ProblemData *PD,
                  const Site a, const Site b, /* by value? */
                  Vec pot)
{
  /* Pair   interaction  parameters:  geometric   average,  arithmetic
     average, and charge product. Site coordinates are not used. */
  const real epsilon = sqrt (a.epsilon * b.epsilon);
  const real sigma = 0.5 * (a.sigma + b.sigma);
  const real q2 = a.charge * b.charge;

  real interval[2];
  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  real h[3];
  FOR_DIM
    h[dim] = PD->h[dim];

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
          real r[3];
          FOR_DIM
            r[dim] = i[dim] * h[dim] + interval[0];

          const real r_s = sqrt (SQR (r[0]) + SQR (r[1]) + SQR (r[2]));

          pot_[i[2]][i[1]][i[0]] +=
            Lennard_Jones (r_s, epsilon, sigma) + Coulomb_short (r_s, q2);
        }
  DAVecRestoreArray (da, pot, &pot_);
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
  There  are several  closure  relations, all  must  satisfy the  same
  interface. Vec c is intent(out) the rest is input.

  typedef void (*Closure)(real beta, Vec v, Vec t, Vec c);

  1)  Hypernetted  Chain  (HNC)  closure relation  to  compute  direct
  correlation function c in real space.  See OZ equation below for the
  second relation between two  unknowns.  The indirect correlation γ =
  h - c is denoted by latin  "t" in other sources. We will use that to
  avoid greek identifiers and confusion with distribution functions:

    c := exp (-βv + γ) - 1 - γ
*/
static void compute_c_HNC (real beta, Vec v, Vec t, Vec c)
{
  real pure f (real v, real t)
  {
    /* exp (-beta * v + t) - 1.0 - t */
    return expm1 (-beta * v + t) - t;
  }
  bgy3d_vec_map2 (c, f, v, t);
}


/* For x <= 0 the same as  exp(x) - 1, but does not grow exponentially
   for positive x:  */
static real lexpm1 (real x)
{
  if (x <= 0.0)
    return expm1 (x);
  else
    return x;
}


/* 2) Kovalenko-Hirata (KH) closure. */
static void compute_c_KH (real beta, Vec v, Vec t, Vec c)
{
  real pure f (real v, real t)
  {
    /* Note that lexpm1() /= expm1(): */
    return lexpm1 (-beta * v + t) - t;
  }
  bgy3d_vec_map2 (c, f, v, t);
}


/*
  3) Percus-Yevick (PY) closure  relation between direct- and indirect
  correlation c and γ:

    c := exp (-βv) [1 + γ] - 1 - γ
*/
static void compute_c_PY (real beta, Vec v, Vec t, Vec c)
{
  real pure f (real v, real t)
  {
    return exp (-beta * v) * (1 + t) - 1 - t;
  }
  bgy3d_vec_map2 (c, f, v, t);
}


static void compute_c (real beta, Vec v, Vec t, Vec c)
{
  char closure[20] = "HNC";
  bgy3d_getopt_string ("--closure", closure, sizeof closure);

  if (strcmp (closure, "HNC") == 0)
    compute_c_HNC (beta, v, t, c);
  else if (strcmp (closure, "KH") == 0)
    compute_c_KH (beta, v, t, c);
  else if (strcmp (closure, "PY") == 0)
    compute_c_PY (beta, v, t, c);
  else
    {
      PetscPrintf (PETSC_COMM_WORLD, "No such OZ closure: %s\n", closure);
      exit (1);
    }
}


static int delta (int i, int j)
{
  return (i == j) ? 1 : 0;
}


/*
  Use the k-representation of Ornstein-Zernike (OZ) equation

    h = c + ρ c * h

  to compute γ =  h - c form c:

                 -1                -1   2
    γ =  (1 - ρc)   c - c = (1 - ρc)  ρc

  If you scale  c by h3 beforehand or  pass rho' = rho *  h3 and scale
  the  result by h3  in addition,  you will  compute exactly  what the
  older version of the function did:
*/
static void compute_t_1 (real rho, Vec c_fft, Vec t_fft)
{
  complex pure f (complex c)
  {
    /* return c / (1.0 - rho * c) - c; */
    return rho * (c * c) / (1.0 - rho * c);
  }
  bgy3d_vec_fft_map1 (t_fft, f, c_fft);
}


/* So far rho is scalar, it could be different for all sites: */
static void compute_t_m (int m, real rho, Vec c_fft[m][m], Vec t_fft[m][m])
{
  complex *c_fft_[m][m], *t_fft_[m][m];

  /* for j <= i only: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* FIXME: vec_get_array() returns real*: */
        c_fft_[i][j] = (void*) vec_get_array (c_fft[i][j]);
        t_fft_[i][j] = (void*) vec_get_array (t_fft[i][j]);

        assert (c_fft[i][j] == c_fft[j][i]);
        assert (t_fft[i][j] == t_fft[j][i]);
      }

  if (m == 0) return;           /* see ref to c_fft[0][0] */

  const int n = vec_local_size (c_fft[0][0]);
  assert (n % 2 == 0);

  {
    complex C[m][m], A[m][m], B[m][m];

    for (int k = 0; k < n / 2; k++)
      {
        /* for j <= i only: */
        for (int i = 0; i < m; i++)
          for (int j = 0; j <= i; j++)
            {
              /*
                A := 1 - ρc
                B := c
              */
              C[i][j] = C[j][i] = c_fft_[i][j][k];
              B[i][j] = B[j][i] = C[i][j];
              A[i][j] = A[j][i] = delta (i, j) - rho * C[i][j];
            }

        /*      -1                -1
          B := A   * B == (1 - ρc)   * c
        */
        hnc3d_sles_zgesv (m, A, B);

        /*                     -1    2
          A := B * C = (1 - ρc)   * c
        */
        hnc3d_sles_zgemm (m, B, C, A); /* FIXME: O(m^3)! */

        /*                  -1     2
          t := ρA = (1 - ρc)   * ρc
        */
        for (int i = 0; i < m; i++)
          for (int j = 0; j <= i; j++)
            t_fft_[i][j][k] = rho * A[i][j];
      }
  }

  /* for j <= i only: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* FIXME: VecRestoreArray() expects real***: */
        VecRestoreArray (c_fft[i][j], (void*) &c_fft_[i][j]);
        VecRestoreArray (t_fft[i][j], (void*) &t_fft_[i][j]);
      }
}


static void compute_t (int m, real rho, Vec c_fft[m][m], Vec t_fft[m][m])
{
  if (m == 1)
    compute_t_1 (rho, c_fft[0][0], t_fft[0][0]); /* faster */
  else
    compute_t_m (m, rho, c_fft, t_fft); /* works for any m */
}


/*
  HNC iteration for an indirect correlation γ = h - c (here denoted by
  t) where other intermediates are considered a function of that:

  t    ->  dt = t    - t
    in           out    in
*/
typedef struct Ctx_t2
{
  State *HD;
  Vec v, c;                     /* real */
  Vec t_fft, c_fft;             /* complex */
} Ctx_t2;


static void iterate_t2 (Ctx_t2 *ctx, Vec t, Vec dt)
{
  const State *HD = ctx->HD;
  const ProblemData *PD = HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real L = PD->interval[1] - PD->interval[0];
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  compute_c (beta, ctx->v, t, ctx->c);

  MatMult (HD->fft_mat, ctx->c, ctx->c_fft);

  /* Translate distribution to the grid corner. */
  bgy3d_vec_fft_trans (HD->dc, PD->N, ctx->c_fft);

  VecScale (ctx->c_fft, h3);

  /* FIXME: compute_t() expects Vec[m][m]: */
  compute_t (1, rho, (void*) &ctx->c_fft, (void*) &ctx->t_fft);

  /* Translate distribution to the grid center. */
  bgy3d_vec_fft_trans (HD->dc, PD->N, ctx->t_fft);

  MatMultTranspose (HD->fft_mat, ctx->t_fft, dt);

  VecScale (dt, 1.0/L/L/L);

  VecAXPY (dt, -1.0, t);
}


/*
  HNC iteration for a direct correlation where the total correlation
  is considered a function of that:

  c    ->  dc = c    - c
    in           out    in
*/
typedef struct Ctx_c2
{
  State *HD;
  int m;
  Vec *v, *t;                   /* [m][m], real */
  Vec *t_fft, *c_fft;           /* [m][m], complex */
} Ctx_c2;


static void iterate_c2 (Ctx_c2 *ctx, Vec C, Vec dC)
{
  /* aliases of the right shape: */
  const int m = ctx->m;
  Vec (*c_fft)[m] = (void*) ctx->c_fft;
  Vec (*t_fft)[m] = (void*) ctx->t_fft;
  Vec (*t)[m] = (void*) ctx->t;
  Vec (*v)[m] = (void*) ctx->v;

  const State *HD = ctx->HD;
  const ProblemData *PD = HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real L = PD->interval[1] - PD->interval[0];
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  /* Establish aliases to the subsections of the long Vec C and dC: */
  local Vec c[m][m], dc[m][m];
  bgy3d_vec_aliases_create2 (C, m, c);
  bgy3d_vec_aliases_create2 (dC, m, dc);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        MatMult (HD->fft_mat, c[i][j], c_fft[i][j]);

        /* Translate distribution to the grid corner. */
        bgy3d_vec_fft_trans (HD->dc, PD->N, c_fft[i][j]);

        VecScale (c_fft[i][j], h3);
      }

  /* Solves the linear equation for t_fft[][]: */
  compute_t (m, rho, c_fft, t_fft);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* Translate distribution to the grid center. */
        bgy3d_vec_fft_trans (HD->dc, PD->N, t_fft[i][j]);

        MatMultTranspose (HD->fft_mat, t_fft[i][j], t[i][j]);

        VecScale (t[i][j], 1.0/L/L/L);

        compute_c (beta, v[i][j], t[i][j], dc[i][j]);
      }

  /* This  destroys the  aliases, but  does  not free  the memory,  of
     course. The actuall data is owned by Vec C and Vec dC: */
  bgy3d_vec_aliases_destroy2 (m, c);
  bgy3d_vec_aliases_destroy2 (m, dc);

  VecAXPY (dC, -1.0, C);
}


/*
  Solving for indirect  correlation γ = h - c  and other quantities of
  HNC equation.  Indirect correlation  γ appears as a primary variable
  here:
*/
static void solvent_solve_t2 (const ProblemData *PD,
                              int m, const Site solvent[m],
                              Vec g[m][m])
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  PetscPrintf (PETSC_COMM_WORLD, "(iterations for γ)\n");

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */
  local Vec c = bgy3d_vec_create (HD->da);
  local Vec c_fft = bgy3d_vec_create (HD->dc); /* complex */
  local Vec t_fft = bgy3d_vec_create (HD->dc); /* complex */

  local Vec v[m][m];
  bgy3d_vec_create2 (HD->da, m, v); /* solvent-solvent interaction */

  /* Get solvent-solvent site-site interactions: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      pair (HD->da, HD->PD, solvent[i], solvent[j], v[i][j]);

  /* Create intial guess: */
  Vec t = bgy3d_vec_create (HD->da); /* FIXME: not deallocated */
  VecSet (t, 0.0);

  assert (m == 1);
  /*
    Find a t  such that dt as returned by iterate_t2  (&ctx, t, dt) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: struct Ctx_t2* vs. void*:
  */
  {
    Ctx_t2 ctx =
      {
        .HD = HD,
        .v = v[0][0],
        .c = c,
        .t_fft = t_fft,
        .c_fft = c_fft,
      };
    bgy3d_snes_default (PD, &ctx, (Function) iterate_t2, t);
  }

  bgy3d_vec_save ("t00.bin", t);

  /* g = γ + c + 1, store in Vec t: */
  VecAXPY (t, 1.0, c);
  VecShift (t, 1.0);

  bgy3d_vec_save ("c00.bin", c);
  bgy3d_vec_save ("g00.bin", t);

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&t); */
  bgy3d_vec_destroy (&c);
  bgy3d_vec_destroy (&c_fft);
  bgy3d_vec_destroy (&t_fft);
  bgy3d_vec_destroy2 (m, v);

  bgy3d_state_destroy (HD);

  /* Return just one distribution. bgy3d_vec_destroy (&) it! */
  g[0][0] = t;
}


/* Solving  for c  and h(c)  of  HNC equation.   Direct correlation  c
   appears as a primary variable here: */
static void solvent_solve_c2 (const ProblemData *PD,
                              int m, const Site solvent[m],
                              Vec g[m][m])
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */

  Vec t[m][m];                  /* FIXME: not deallocated */
  bgy3d_vec_create2 (HD->da, m, t);

  local Vec t_fft[m][m];
  bgy3d_vec_create2 (HD->dc, m, t_fft); /* complex */

  local Vec c_fft[m][m];
  bgy3d_vec_create2 (HD->dc, m, c_fft); /* complex */

  local Vec v[m][m];
  bgy3d_vec_create2 (HD->da, m, v); /* solvent-solvent interaction */

  /* Get solvent-solvent site-site interactions: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      pair (HD->da, HD->PD, solvent[i], solvent[j], v[i][j]);

  /*
    For primary  variable there are two  ways to access  the data: via
    the long  Vec and  m * (m  + 1)  / 2 shorter  Vecs aliased  to the
    subsections of the longer one.
  */
  local Vec C = bgy3d_vec_pack_create2 (HD->da, m); /* long Vec */

  local Vec c[m][m];
  bgy3d_vec_aliases_create2 (C, m, c); /* aliases to subsections */

  /* Create intial guess: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecSet (c[i][j], 0.0);

  /*
    Find a C  such that dC as returned by iterate_c2  (&ctx, C, dC) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: struct Ctx_c2* vs. void*:
  */
  {
    Ctx_c2 ctx =
      {
        .HD = HD,
        .m = m,
        .v = (void*) v,
        .t = (void*) t,
        .t_fft = (void*) t_fft,
        .c_fft = (void*) c_fft,
      };
    bgy3d_snes_default (PD, &ctx, (Function) iterate_c2, C);
  }

  bgy3d_vec_save2 ("t%d%d.bin", m, t);

  /* g = γ + c + 1, store in Vec t: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        VecAXPY (t[i][j], 1.0, c[i][j]);
        VecShift (t[i][j], 1.0);
      }

  bgy3d_vec_save2 ("c%d%d.bin", m, c);
  bgy3d_vec_save2 ("g%d%d.bin", m, t);

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&t); */
  bgy3d_vec_destroy2 (m, t_fft);
  bgy3d_vec_destroy2 (m, c_fft);
  bgy3d_vec_destroy2 (m, v);

  bgy3d_vec_aliases_destroy2 (m, c);
  bgy3d_vec_pack_destroy2 (&C);

  bgy3d_state_destroy (HD);

  /* Dont forget to bgy3d_vec_destroy2() it! */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      g[i][j] = g[j][i] = t[i][j];
}


/* Solving  for c  and h(c)  of  HNC equation.   Direct correlation  c
   appears as a primary variable here: */
void hnc3d_solvent_solve (const ProblemData *PD,
                          int m, const Site solvent[m],
                          Vec g[m][m])
{
  if (1)
    solvent_solve_c2 (PD, m, solvent, g);
  else
    solvent_solve_t2 (PD, m, solvent, g);
}


/*
  Solving  for c  and  h(c)  of HNC  equation  with non-linear  solver
  specified by  the --snes-solver flag.  Direct  correlation c appears
  as a primary variable here:
*/
Vec HNC3d_solvent_solve (const ProblemData *PD, Vec g_ini)
{
  assert (g_ini == PETSC_NULL);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving 3d-HNC equation.\n");

  int m;                        /* number of solvent sites */
  const Site *solvent;          /* solvent[m] */

  char name[200] = "LJ";        /* solvent & default solute */

  /* Get the number of solvent sites and their parameters. Get it from
     the solute tables: */
  bgy3d_solute_get (name, &m, &solvent);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", name);

  Vec g[m][m];
  hnc3d_solvent_solve (PD, m, solvent, g);

  assert (m == 1);
  return g[0][0];               /* bgy3d_vec_destroy (&) it! */
}


/*
  HNC iteration for a fixed direct correlation:

  h    ->  dh = h    - h
    in           out    in
*/
typedef struct Ctx_h1
{
  State *HD;
  Vec v, t;                     /* real */
  Vec c_fft, h_fft, t_fft;      /* complex */
} Ctx_h1;


static void iterate_h1 (Ctx_h1 *ctx, Vec h, Vec dh)
{
  const ProblemData *PD = ctx->HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real N3 = PD->N[0] * PD->N[1] * PD->N[2];

  /* fft(h).  Here h is the 3d  unknown hole density h1 of the solvent
     sites. */
  MatMult (ctx->HD->fft_mat, h, ctx->h_fft);

  /* fft(h)  *  fft(c).  Here   c  is  the  constant  (radial)  direct
     correlation c2 of the pure solvent. */
  complex pure mul (complex x, complex y)
  {
    return x * y;
  }
  bgy3d_vec_fft_map2 (ctx->t_fft, mul, ctx->c_fft, ctx->h_fft);

  /* t = fft^-1 (fft(c) * fft(h)). Here t is 3d t1. */
  MatMultTranspose (ctx->HD->fft_mat, ctx->t_fft, ctx->t);

  VecScale (ctx->t, rho * h3 / N3);

  /*
    The new candidate for the total correlation

      h = exp (-βv + t) - 1

    with

      t = ρ (c * h)

    computed using the input h.  Resulting new h will be stored in Vec
    dh. To compute that we re-use the closure relation, which, in case
    of HNC closure, outputs

      exp (-βv + t) - 1.0 - t

    (I dont call the expression above "c" as in the closure definition
    to avoid the necessity to differentiate between pair and singleton
    c2 and c1).   It is questionable whether we  should always use HNC
    closure  here  or  respect   the  input  of  the  user  (--closure
    switch). There is  only a chance for solute  == solvent to produce
    the same result as the pure solvent if we use the "native" closure
    here:
  */
  compute_c (beta, ctx->v, ctx->t, dh);
  VecAXPY (dh, 1.0, ctx->t);

  /*
    dh := h    - h
           out    in
  */
  VecAXPY (dh, -1.0, h);
}


/*
  HNC iteration for a fixed direct correlation:

  t    ->  dt = t    - t
    in           out    in
*/
typedef struct Ctx_t1
{
  State *HD;
  Vec v, h;                     /* real */
  Vec c_fft, h_fft, t_fft;      /* complex */
} Ctx_t1;


static void iterate_t1 (Ctx_t1 *ctx, Vec t, Vec dt)
{
  const ProblemData *PD = ctx->HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real N3 = PD->N[0] * PD->N[1] * PD->N[2];

  /*
    The new candidate for the total correlation

      h = exp (-βv + t) - 1

    computed using the input t.  To compute that we re-use the closure
    relation, which, in case of HNC closure, outputs

      exp (-βv + t) - 1.0 - t

    See comments on the use of closure in iterate_h1().
  */
  compute_c (beta, ctx->v, t, ctx->h);
  VecAXPY (ctx->h, 1.0, t);

  /* fft(h).  Here h is the 3d  unknown hole density h1 of the solvent
     sites. */
  MatMult (ctx->HD->fft_mat, ctx->h, ctx->h_fft);

  /* fft(h)  *  fft(c).  Here   c  is  the  constant  (radial)  direct
     correlation c2 of the pure solvent. */
  complex pure mul (complex x, complex y)
  {
    return x * y;
  }
  bgy3d_vec_fft_map2 (ctx->t_fft, mul, ctx->c_fft, ctx->h_fft);

  /* t = fft^-1 (fft(c) * fft(h)). Here t is 3d t1. */
  MatMultTranspose (ctx->HD->fft_mat, ctx->t_fft, dt);

  VecScale (dt, rho * h3 / N3);

  /*
    dt := t    - t
           out    in
  */
  VecAXPY (dt, -1.0, t);
}


/* Reads c2[][] as previousely written by solvent solver. */
static void solvent_kernel (State *HD, int m, Vec c_fft[m][m])
{
  local Vec c[m][m];
  bgy3d_vec_create2 (HD->da, m, c);

  /* Load c_1d from file: */
  if (bgy3d_getopt_test ("--from-radial-g2")) /* FIXME: better name? */
    bgy3d_vec_read_radial2 (HD->da, HD->PD, "c%d%d.txt", m, c);
  else
    bgy3d_vec_read2 ("c%d%d.bin", m, c);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        MatMult (HD->fft_mat, c[i][j], c_fft[i][j]);

        /*
          Translate the distribution to  the grid corner. This is what
          one expects  in convolution  integrals. FIXME: or  should we
          rather store the convolution kernel on disk in ready form?
        */
        bgy3d_vec_fft_trans (HD->dc, HD->PD->N, c_fft[i][j]);
      }
  bgy3d_vec_destroy2 (m, c);
}


/* Solving for h only of HNC equation. Direct correlation c appears as
   an input here */
static void solute_solve_h1 (const ProblemData *PD,
                             const int m, const Site solvent[m],
                             const int n, const Site solute[n],
                             Vec g[m])
{
  assert (m == 1);

  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */
  local Vec t = bgy3d_vec_create (HD->da);
  local Vec c_fft = bgy3d_vec_create (HD->dc);  /* complex */
  local Vec h_fft = bgy3d_vec_create (HD->dc);  /* complex */
  local Vec t_fft = bgy3d_vec_create (HD->dc);  /* complex */

  local Vec v[m];
  bgy3d_vec_create1 (HD->da, m, v); /* solute-solvent interaction */

  /*
    Get  the solvent-solvent direct  correlation function.   FIXME: we
    get the solvent  description as an argument, but  read the file to
    get  the  rest.   There  is  no  guarantee  the  two  sources  are
    consistent.
  */
  solvent_kernel (HD, 1, (void*) &c_fft);

  /*
    Get  solute-solvent interaction.   Fill v[]  with  the short-range
    potential. FIXME: long-range Coulomb is ignored yet:
  */
  bgy3d_solute_field (HD, m, solvent, n, solute,
                      v,          /* out */
                      NULL, NULL, /* no coulomb */
                      NULL);      /* no electrons */

  assert (m == 1);

  /* Create global vectors */
  Vec h = bgy3d_vec_duplicate (v[0]); /* FIXME: not deallocated */

  /* Set initial guess */
  compute_h (HD->PD->beta, v[0], h);

  /*
    Find  an h such  that dh  as returned  by iterate  (HD, h,  dh) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: Ctx_h1* vs. void*:
  */
  {
    /* Work area for iterate_h1(): */
    Ctx_h1 ctx =
      {
        .HD = HD,
        .v = v[0],
        .t = t,
        .c_fft = c_fft,
        .h_fft = h_fft,
        .t_fft = t_fft,
      };

    bgy3d_snes_default (PD, &ctx, (Function) iterate_h1, h);
  }

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&h) */
  bgy3d_vec_destroy (&t);
  bgy3d_vec_destroy (&h_fft);
  bgy3d_vec_destroy (&c_fft);
  bgy3d_vec_destroy (&t_fft);
  bgy3d_vec_destroy1 (m, v);

  bgy3d_state_destroy (HD);

  /* g := h + 1 */
  VecShift (h, 1.0);

  bgy3d_vec_save ("g0.bin", h);

  /* Just one distribution so far. bgy3d_vec_destroy (&) it! */
  g[0] = h;
}


/* Solving for h only of HNC equation. Direct correlation c appears as
   an input here */
static void solute_solve_t1 (const ProblemData *PD,
                             const int m, const Site solvent[m],
                             const int n, const Site solute[n],
                             Vec g[m])
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */

  local Vec t[m];
  bgy3d_vec_create1 (HD->da, m, t);

  local Vec h_fft[m];
  bgy3d_vec_create1 (HD->dc, m, h_fft); /* complex */

  local Vec t_fft[m];
  bgy3d_vec_create1 (HD->dc, m, t_fft); /* complex */

  local Vec v[m];
  bgy3d_vec_create1 (HD->da, m, v); /* solute-solvent interaction */

  /* This should be the only pair quantity: */
  local Vec c_fft[m][m];
  bgy3d_vec_create2 (HD->dc, m, c_fft);        /* complex */

  /*
    Get  the solvent-solvent direct  correlation function.   FIXME: we
    get the solvent  description as an argument, but  read the file to
    get  the  rest.   There  is  no  guarantee  the  two  sources  are
    consistent.
  */
  solvent_kernel (HD, m, c_fft);

  /*
    Get  solute-solvent interaction.   Fill v[]  with  the short-range
    potential. FIXME: long-range Coulomb is ignored yet:
  */
  bgy3d_solute_field (HD, m, solvent, n, solute,
                      v,          /* out */
                      NULL, NULL, /* no coulomb */
                      NULL);      /* no electrons */

  /* Create global vectors */
  Vec h[m];
  bgy3d_vec_create1 (HD->da, m, h); /* FIXME: not deallocated */

  /* Set initial guess */
  for (int i = 0; i < m; i++)
    compute_h (HD->PD->beta, v[i], h[i]);

  assert (m == 1);

  /*
    Find a  t such that  dt as returned  by iterate_t1 (HD, t,  dt) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: Ctx_t1* vs. void*:
  */
  {
    /* Work area for iterate_h1(): */
    Ctx_t1 ctx =
      {
        .HD = HD,
        .v = v[0],
        .h = h[0],
        .c_fft = c_fft[0][0],
        .h_fft = h_fft[0],
        .t_fft = t_fft[0],
      };

    bgy3d_snes_default (PD, &ctx, (Function) iterate_t1, t[0]);
  }
  /*
    FIXME: by ignoring the value of  Vec t and using Vec h instead, we
    assume  that the two  are consistent.   Hopefully the  SNES solver
    does not  return "untested" solutions.  A call  to iterate_t1() to
    "test" a solution would update Vec h, namely.
  */

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&h) */
  bgy3d_vec_destroy1 (m, t);
  bgy3d_vec_destroy1 (m, h_fft);
  bgy3d_vec_destroy1 (m, t_fft);
  bgy3d_vec_destroy1 (m, v);

  /* This should be the only pair quantity: */
  bgy3d_vec_destroy2 (m, c_fft);

  bgy3d_state_destroy (HD);

  /* g := h + 1 */
  for (int i = 0; i < m; i++)
    VecShift (h[i], 1.0);

  bgy3d_vec_save1 ("g%d.bin", m, h);

  /* Caller is supposed to destroy it! */
  for (int i = 0; i < m; i++)
    g[i] = h[i];
}


void hnc3d_solute_solve (const ProblemData *PD,
                         const int m, const Site solvent[m],
                         const int n, const Site solute[n],
                         Vec g[m])
{
  /*
    At least for water in OW solvent (water-like model, just a neutral
    oxygen with LJ parameters of  water) Newton iteration over t = log
    [g exp(βv)] converges, whereas iteration over  h = g - 1 does not.
    The problem with  h as a primary variable cannot  be severe as the
    Jager solver manages that.
  */
  if (0)
    solute_solve_h1 (PD, m, solvent, n, solute, g);
  else
    solute_solve_t1 (PD, m, solvent, n, solute, g);
}


/*
  Solving for  h of HNC equation  with a default  non-linear solver as
  specified  by  the  --snes-solver  option.  The  direct  correlation
  function "c" is fixed and appears as an input here:
*/
Vec HNC3d_solute_solve (const ProblemData *PD, Vec g_ini)
{
  assert (g_ini == PETSC_NULL);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving 3d-HNC equation. Fixed c.\n");

  int n;                        /* number of solute sites */
  const Site *solute;           /* solute[n] */

  int m;                        /* number of solvent sites */
  const Site *solvent;          /* solvent[m] */

  char name[200] = "LJ";        /* solvent & default solute */

  /* Get the number of solvent sites and their parameters. Get it from
     the solute tables: */
  bgy3d_solute_get (name, &m, &solvent);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", name);

  /* Get solute name or stay with the default: */
  bgy3d_getopt_string ("--solute", name, sizeof name);

  /* Get the solute from the tables: */
  bgy3d_solute_get (name, &n, &solute);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solute is %s.\n", name);

  Vec g[m];
  hnc3d_solute_solve (PD, m, solvent, n, solute, g);

  assert (m == 1);
  return g[0];
}
