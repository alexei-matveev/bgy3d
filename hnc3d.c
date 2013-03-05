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


void compute_c (real beta, Vec v, Vec t, Vec c)
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


#ifdef HNC3D_T
/*
  HNC iteration for an indirect correlation γ = h - c (here denoted by
  t) where other intermediates are considered a function of that:

  t    ->  dt = t    - t
    in           out    in
*/
typedef struct Ctx_t
{
  State *HD;
  Vec v, c;                     /* real */
  Vec t_fft, c_fft;             /* complex */
} Ctx_t;


static void iterate_t (Ctx_t *ctx, Vec t, Vec dt)
{
  const State *HD = ctx->HD;
  const ProblemData *PD = HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real L = PD->interval[1] - PD->interval[0];
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  compute_c (beta, ctx->v, t, ctx->c);

  MatMult (HD->fft_mat, ctx->c, ctx->c_fft);

  VecScale (ctx->c_fft, h3);

  compute_t (rho, ctx->c_fft, ctx->t_fft);

  MatMultTranspose (HD->fft_mat, ctx->t_fft, dt);

  VecScale (dt, 1.0/L/L/L);

  VecAXPY (dt, -1.0, t);
}

#else
/*
  HNC iteration for a direct correlation where the total correlation
  is considered a function of that:

  c    ->  dc = c    - c
    in           out    in
*/
typedef struct Ctx_c
{
  State *HD;
  Vec v, t;                     /* real */
  Vec t_fft, c_fft;             /* complex */
} Ctx_c;


static void iterate_c (Ctx_c *ctx, Vec c, Vec dc)
{
  const State *HD = ctx->HD;
  const ProblemData *PD = HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real L = PD->interval[1] - PD->interval[0];
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  MatMult (HD->fft_mat, c, ctx->c_fft);

  VecScale (ctx->c_fft, h3);

  compute_t (rho, ctx->c_fft, ctx->t_fft);

  MatMultTranspose (HD->fft_mat, ctx->t_fft, ctx->t);

  VecScale (ctx->t, 1.0/L/L/L);

  compute_c (beta, ctx->v, ctx->t, dc);

  VecAXPY (dc, -1.0, c);
}
#endif

#ifdef HNC3D_T
/*
  Solving for indirect  correlation γ = h - c  and other quantities of
  HNC equation.  Indirect correlation  γ appears as a primary variable
  here:
*/
void hnc3d_solvent_solve (const ProblemData *PD, Vec g[1][1])
{
  PetscPrintf (PETSC_COMM_WORLD, "(iterations for γ)\n");

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */
  Vec v = bgy3d_vec_create (HD->da); /* solvent-solvent interaction */
  Vec c = bgy3d_vec_create (HD->da);
  Vec c_fft = bgy3d_vec_create (HD->dc); /* complex */
  Vec t_fft = bgy3d_vec_create (HD->dc); /* complex */


  /* FIXME:   this  is   abused  to   get  both   solvent-solvent  and
     solute-solvent interactions: */
  solute_field (HD->da, HD->PD, v);

  /* Create intial guess: */
  Vec t = bgy3d_vec_create (HD->da);
  VecSet (t, 0.0);

  /*
    Find a  t such that dt as  returned by iterate_t (&ctx,  t, dt) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: struct Ctx_t* vs. void*:
  */
  {
    Ctx_t ctx =
      {
        .HD = HD,
        .v = v,
        .c = c,
        .t_fft = t_fft,
        .c_fft = c_fft,
      };
    bgy3d_snes_default (PD, &ctx, (Function) iterate_t, t);
  }

  bgy3d_vec_save ("t00.bin", t);

  /* g = γ + c + 1, store in Vec t: */
  VecAXPY (t, 1.0, c);
  VecShift (t, 1.0);

  bgy3d_vec_save ("c00.bin", c);
  bgy3d_vec_save ("g00.bin", t);

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&t); */
  bgy3d_vec_destroy (&v);
  bgy3d_vec_destroy (&c);
  bgy3d_vec_destroy (&c_fft);
  bgy3d_vec_destroy (&t_fft);

  bgy3d_state_destroy (HD);

  /* Return just one distribution. bgy3d_vec_destroy (&) it! */
  g[0][0] = t;
}

#else
/* Solving  for c  and h(c)  of  HNC equation.   Direct correlation  c
   appears as a primary variable here: */
void hnc3d_solvent_solve (const ProblemData *PD,
                          int m, const Site solvent[m],
                          Vec g[m][m])
{
  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */
  Vec t = bgy3d_vec_create (HD->da);
  Vec t_fft = bgy3d_vec_create (HD->dc); /* complex */
  Vec c_fft = bgy3d_vec_create (HD->dc); /* complex */

  Vec v[m][m];
  bgy3d_vec_create2 (HD->da, m, v); /* solvent-solvent interaction */

  /* Get solvent-solvent site-site interactions: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      pair (HD->da, HD->PD, solvent[i], solvent[j], v[i][j]);

  /* Create intial guess: */
  Vec c = bgy3d_vec_create (HD->da);
  VecSet (c, 0.0);

  assert (m == 1);
  /*
    Find a  c such that dc as  returned by iterate_c (&ctx,  c, dc) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: struct Ctx_c* vs. void*:
  */
  {
    Ctx_c ctx =
      {
        .HD = HD,
        .v = v[0][0],
        .t = t,
        .t_fft = t_fft,
        .c_fft = c_fft,
      };
    bgy3d_snes_default (PD, &ctx, (Function) iterate_c, c);
  }

  /* g = γ + c + 1, store in Vec t: */
  VecAXPY (t, 1.0, c);
  VecShift (t, 1.0);

  bgy3d_vec_save ("c00.bin", c);
  bgy3d_vec_save ("g00.bin", t);

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&t); */
  bgy3d_vec_destroy (&t_fft);
  bgy3d_vec_destroy (&c);
  bgy3d_vec_destroy (&c_fft);
  bgy3d_vec_destroy2 (m, v);

  bgy3d_state_destroy (HD);

  /* Return just one distribution. bgy3d_vec_destroy (&) it! */
  g[0][0] = t;
}
#endif

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
typedef struct Ctx_h
{
  State *HD;
  Vec v, t;                     /* real */
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

    h = exp (-βv + ρt) - 1

    stored in Vec dh:
  */
  real pure h_out (real v, real t)
  {
    return expm1 (-beta * v + rho * t);
  }
  bgy3d_vec_map2 (dh, h_out, ctx->v, t);

  /*
    dh := h    - h
           out    in
  */
  VecAXPY (dh, -1.0, h);
}


static void solvent_kernel (State *HD, Vec c_fft)
{
  Vec c = bgy3d_vec_create (HD->da);

  /* Load c_1d from file: */
  if (bgy3d_getopt_test ("--from-radial-g2")) /* FIXME: better name? */
    bgy3d_vec_read_radial (HD->da, HD->PD, "c1dfile", c);
  else
    bgy3d_vec_read ("c00.bin", c);

  MatMult (HD->fft_mat, c, c_fft);

  bgy3d_vec_destroy (&c);
}


/* Solving for h only of HNC equation. Direct correlation c appears as
   an input here */
void hnc3d_solute_solve (const ProblemData *PD,
                         const int m, const Site solvent[m],
                         const int n, const Site solute[n],
                         Vec g[m])
{
  assert (m == 1);

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */
  Vec t = bgy3d_vec_create (HD->da);
  Vec c_fft = bgy3d_vec_create (HD->dc);  /* complex */
  Vec h_fft = bgy3d_vec_create (HD->dc);  /* complex */
  Vec ch_fft = bgy3d_vec_create (HD->dc); /* complex */

  Vec v[m];
  bgy3d_vec_create1 (HD->da, m, v); /* solute-solvent interaction */

  /*
    Get  the solvent-solvent direct  correlation function.   FIXME: we
    get the solvent  description as an argument, but  read the file to
    get  the  rest.   There  is  no  guarantee  the  two  sources  are
    consistent.
  */
  solvent_kernel (HD, c_fft);

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
  Vec h = bgy3d_vec_duplicate (v[0]);

  /* Set initial guess */
  compute_h (HD->PD->beta, v[0], h);

  /*
    Find  an h such  that dh  as returned  by iterate  (HD, h,  dh) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: Ctx_h* vs. void*:
  */
  {
    /* Work area for iterate_h(): */
    Ctx_h ctx =
      {
        .HD = HD,
        .v = v[0],
        .t = t,
        .c_fft = c_fft,
        .h_fft = h_fft,
        .ch_fft = ch_fft,
      };

    bgy3d_snes_default (PD, &ctx, (Function) iterate_h, h);
  }

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&h) */
  bgy3d_vec_destroy (&t);
  bgy3d_vec_destroy (&h_fft);
  bgy3d_vec_destroy (&c_fft);
  bgy3d_vec_destroy (&ch_fft);
  bgy3d_vec_destroy1 (m, v);

  bgy3d_state_destroy (HD);

  /* g := h + 1 */
  VecShift (h, 1.0);

  bgy3d_vec_save ("g0.bin", h);

  /* Just one distribution so far. bgy3d_vec_destroy (&) it! */
  g[0] = h;
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
