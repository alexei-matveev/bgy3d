/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-vec.h"          /* bgy3d_vec_create() */
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-force.h"        /* bgy3d_pair_potential() */
#include "bgy3d-pure.h"         /* bgy3d_omega_fft_create() */
#include "bgy3d-snes.h"         /* bgy3d_snes_default() */
#include "hnc3d-sles.h"         /* hnc3d_sles_zgesv() */
#include "hnc3d.h"
#include <math.h>               /* expm1() */


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
  vec_map2 (c, f, v, t);
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
  vec_map2 (c, f, v, t);
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
  vec_map2 (c, f, v, t);
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


/*
  Compute

    h = c'(t) + t

  In particular  for t =  0 returns  h = exp(-βv)  - 1 (except  for KH
  closure).   Used  to  compute   the  new  candidate  for  the  total
  correlation with the original expression being

    h = exp (-βv + t) - 1

  where

    t = ρ (c *  h)

  is computed  using the current input  h in HNC3d  model.  To compute
  that we re-use the closure  relation, which, in case of HNC closure,
  outputs

    exp (-βv + t) - 1.0 - t

  (I am somewhat reluctant to call  the expression above "c" as in the
  closure definition  to avoid the necessity  to differentiate between
  pair and singleton c = c2  and c' = c1).  It is questionable whether
  we should  always use HNC closure  here or respect the  input of the
  user  (--closure switch).  There  is  only a  chance  for solute  ==
  solvent to produce the same result as the pure solvent if we use the
  "native" closure here:
*/
static void compute_h (real beta, Vec v, Vec t, Vec h)
{
  compute_c (beta, v, t, h);
  VecAXPY (h, 1.0, t);
}


/* h := exp (-β v) - 1. A kind of special case of the above. */
static void compute_h0 (real beta, Vec v, Vec h)
{
  real pure f (real v)
  {
    /* exp (-beta * v) - 1 */
    return expm1 (-beta * v);
  }
  vec_map1 (h, f, v);
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
static void compute_t2_1 (real rho, Vec c_fft, Vec t_fft)
{
  complex pure f (complex c)
  {
    /* return c / (1.0 - rho * c) - c; */
    return rho * (c * c) / (1.0 - rho * c);
  }
  vec_fft_map1 (t_fft, f, c_fft);
}


/* So far rho is scalar, it could be different for all sites: */
static void
compute_t2_m (int m, real rho, Vec c_fft[m][m], Vec w_fft[m][m], Vec t_fft[m][m])
{
  if (m == 0) return;           /* see ref to c_fft[0][0] */

  local complex *c_fft_[m][m], *w_fft_[m][m], *t_fft_[m][m];

  /* for j <= i only: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* FIXME: vec_get_array() returns real*: */
        c_fft_[i][j] = (void*) vec_get_array (c_fft[i][j]);
        t_fft_[i][j] = (void*) vec_get_array (t_fft[i][j]);

        /* Diagonals are NULL: */
        w_fft_[i][j] = (i == j) ? NULL : (void*) vec_get_array (w_fft[i][j]);

        assert (c_fft[i][j] == c_fft[j][i]);
        assert (t_fft[i][j] == t_fft[j][i]);
        assert (w_fft[i][j] == w_fft[j][i]);
      }

  const int n = vec_local_size (c_fft[0][0]);
  assert (n % 2 == 0);

  {
    complex H[m][m], C[m][m], W[m][m], WC[m][m], T[m][m];

    for (int k = 0; k < n / 2; k++)
      {
        /* for j <= i only: */
        for (int i = 0; i < m; i++)
          for (int j = 0; j <= i; j++)
            {
              /*
                Extract C  and W for  this particular momentum  k from
                scattered arrays into contiguous matrices.
              */
              C[i][j] = C[j][i] = c_fft_[i][j][k];

              /* Diagonal is implicitly 1: */
              W[i][j] = W[j][i] = (i == j) ? 1.0 : w_fft_[i][j][k];
            }

        /* WC,  an  intermediate.  See  comment  on  layout of  result
           below.  The input is symmetric, though: */
        hnc3d_sles_zgemm (m, W, C, WC);

        /* T :=  1 - ρWC.  The  matrix T is  used here as a  free work
           array: */
        for (int i = 0; i < m; i++)
          for (int j = 0; j < m; j++)
            T[i][j] = delta (i, j) - rho * WC[i][j];

        /*
          H := WCW.  Temporarily ---  it will be overwritten with the
          real H  after solving the  linear equations.

          Note that  the matrix multiplication is  performed using the
          Fortran  (column major) interpetation  of the  matrix memory
          layout. In Fortran notation:

            H(i, j) = Σ  WC(i, k) * W(k, j)
                       k

          or, in C-notation:

            H[j][i] = Σ  WC[k][i] * W[j][k]
                       k
        */
        hnc3d_sles_zgemm (m, WC, W, H);

        /*
          Solving the linear equation makes H have the literal meaning
          of the total correlation matrix (input T is destroyed):

                -1                -1
          H := T   * H == (1 - ρc)   * c
        */
        hnc3d_sles_zgesv (m, T, H);


        /*
          The  same  effect  is  achieved  in  1x1  version  of  the  code
          differently:

          T := H - C
        */
        for (int i = 0; i < m; i++)
          for (int j = 0; j <= i; j++)
            t_fft_[i][j][k] = H[i][j] - C[i][j];
      }
  }

  /* for j <= i only: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* FIXME:   vec_restore_array()  expects   real**,   we  offer
           complex** instead: */
        vec_restore_array (c_fft[i][j], (void*) &c_fft_[i][j]);
        vec_restore_array (t_fft[i][j], (void*) &t_fft_[i][j]);
        if (i != j)
          vec_restore_array (w_fft[i][j], (void*) &w_fft_[i][j]);
      }
}


static void
compute_t2 (int m, real rho, Vec c_fft[m][m], Vec w_fft[m][m], Vec t_fft[m][m])
{
  if (m == 1)
    /* Here w = 1, identically: */
    compute_t2_1 (rho, c_fft[0][0], t_fft[0][0]); /* faster */
  else
    compute_t2_m (m, rho, c_fft, w_fft, t_fft); /* works for any m */
}


/*
  There are currently two types of HNC iterations for pure solvent.

  a) HNC iteration for an indirect correlation γ = h - c (here denoted
     by t)

     t    ->  dt = t    - t
      in           out    in

     where other  intermediates are considered a function  of that. In
     particular  in this case  another important  intermediate is  y =
     c(t).

  b) HNC iteration for a direct correlation where the total correlation
     is considered a function of that:

     c    ->  dc = c    - c
      in           out    in

     In this case the indirect  correlation γ is considered a function
     of c and is stored in y = γ(c).
*/
typedef struct Ctx2
{
  State *HD;
  int m;
  Vec *v_short;              /* [m][m], real, center, const */
  Vec *v_long_fft;           /* [m][m], complex, center, const */
  Vec *y;                    /* [m][m], real, either y = c or y = t */
  Vec *t_fft, *c_fft;        /* [m][m], complex */
  Vec *w_fft;                /* [m][m], complex, NULL diagonal */
} Ctx2;


static void iterate_t2 (Ctx2 *ctx, Vec T, Vec dT)
{
  const int m = ctx->m;

  /* aliases of the correct [m][m] shape: */
  Vec (*c_fft)[m] = (void*) ctx->c_fft;
  Vec (*t_fft)[m] = (void*) ctx->t_fft;
  Vec (*w_fft)[m] = (void*) ctx->w_fft;
  Vec (*c)[m] = (void*) ctx->y; /* y = c(t) here */
  Vec (*v_short)[m] = (void*) ctx->v_short;
  Vec (*v_long_fft)[m] = (void*) ctx->v_long_fft;

  const State *HD = ctx->HD;
  const ProblemData *PD = HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real L = PD->interval[1] - PD->interval[0];
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  /* Establish aliases to the subsections of the long Vec T and dT: */
  local Vec t[m][m], dt[m][m];
  bgy3d_vec_aliases_create2 (T, m, t);
  bgy3d_vec_aliases_create2 (dT, m, dt);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        compute_c (beta, v_short[i][j], t[i][j], c[i][j]);

        MatMult (HD->fft_mat, c[i][j], c_fft[i][j]);

        VecScale (c_fft[i][j], h3);

        /*
          The real-space  representation encodes only  the short-range
          part  of  the  direct  corrlation. The  (fixed)  long  range
          contribution is added here:

            C := C  - βV
                  S     L
        */
        VecAXPY (c_fft[i][j], -beta, v_long_fft[i][j]);

        /* Translate distribution to the grid corner. */
        bgy3d_vec_fft_trans (HD->dc, PD->N, c_fft[i][j]);
      }

  /* Solves the OZ linear equation for t_fft[][].: */
  compute_t2 (m, rho, c_fft, w_fft, t_fft);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* Translate distribution to the grid center. */
        bgy3d_vec_fft_trans (HD->dc, PD->N, t_fft[i][j]);

        /*
          Since we plugged in the Fourier transform of the full direct
          correlation  including  the  long  range part  into  the  OZ
          equation what  we get out  is the full  indirect correlation
          including the  long-range part.   The menmonic is  C +  T is
          short range.  Take it out:

            T := T - βV
             S         L
        */
        VecAXPY (t_fft[i][j], -beta, v_long_fft[i][j]);

        MatMultTranspose (HD->fft_mat, t_fft[i][j], dt[i][j]);

        VecScale (dt[i][j], 1.0/L/L/L);
      }

  /*
    This  destroys the  aliases,  but  does not  free  the memory,  of
    course. The actuall data is owned by Vec T and Vec dT. From now on
    one may access T and dT directly again.
  */
  bgy3d_vec_aliases_destroy2 (T, m, t);
  bgy3d_vec_aliases_destroy2 (dT, m, dt);

  VecAXPY (dT, -1.0, T);
}


static void iterate_c2 (Ctx2 *ctx, Vec C, Vec dC)
{
  const int m = ctx->m;

  /* aliases of the correct [m][m] shape: */
  Vec (*c_fft)[m] = (void*) ctx->c_fft;
  Vec (*t_fft)[m] = (void*) ctx->t_fft;
  Vec (*w_fft)[m] = (void*) ctx->w_fft;
  Vec (*t)[m] = (void*) ctx->y; /* y = t(c) here */
  Vec (*v_short)[m] = (void*) ctx->v_short;
  Vec (*v_long_fft)[m] = (void*) ctx->v_long_fft;

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

        VecScale (c_fft[i][j], h3);

        /*
          The real-space  representation encodes only  the short-range
          part  of  the  direct  corrlation. The  (fixed)  long  range
          contribution is added here:

            C := C  - βV
                  S     L
        */
        VecAXPY (c_fft[i][j], -beta, v_long_fft[i][j]);

        /* Translate distribution to the grid corner. */
        bgy3d_vec_fft_trans (HD->dc, PD->N, c_fft[i][j]);
      }

  /* Solves the OZ linear equation for t_fft[][].: */
  compute_t2 (m, rho, c_fft, w_fft, t_fft);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* Translate distribution to the grid center. */
        bgy3d_vec_fft_trans (HD->dc, PD->N, t_fft[i][j]);

        /*
          Since we plugged in the Fourier transform of the full direct
          correlation  including  the  long  range part  into  the  OZ
          equation what  we get out  is the full  indirect correlation
          including the  long-range part.   The menmonic is  C +  T is
          short range.  Take it out:

            T := T - βV
             S         L
        */
        VecAXPY (t_fft[i][j], -beta, v_long_fft[i][j]);

        MatMultTranspose (HD->fft_mat, t_fft[i][j], t[i][j]);

        VecScale (t[i][j], 1.0/L/L/L);

        compute_c (beta, v_short[i][j], t[i][j], dc[i][j]);
      }

  /*
    This  destroys the  aliases,  but  does not  free  the memory,  of
    course. The actuall data is owned by Vec C and Vec dC. From now on
    one may access C and dC directly again.
  */
  bgy3d_vec_aliases_destroy2 (C, m, c);
  bgy3d_vec_aliases_destroy2 (dC, m, dc);

  VecAXPY (dC, -1.0, C);
}


/*
  Solving  for  either  indirect correlation  γ  =  h  - c  or  direct
  correlation c  and other quantities of HNC  equation.  Either direct
  correlation  c  or  indirect  correlation  γ appears  as  a  primary
  variable x here. All other unknowns are functionals of that.
*/
void hnc3d_solvent_solve (const ProblemData *PD,
                          int m, const Site solvent[m],
                          Vec g[m][m])
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */

  /* Flip this to switch between c and t as primary variables: true =>
     c, false => t: */
  const bool yes = false;

  if (yes)
    PetscPrintf (PETSC_COMM_WORLD, "(iterations for c)\n");
  else
    PetscPrintf (PETSC_COMM_WORLD, "(iterations for γ)\n");

  /*
    This will be a functional y(x) of primary variable, either y(x) ==
    t(c)  or  vice  versa, y(x)  =  c(t).  Depending  on the  type  of
    iteration.
  */
  Vec y[m][m];                  /* FIXME: not deallocated */
  bgy3d_vec_create2 (HD->da, m, y);

  local Vec t_fft[m][m];
  bgy3d_vec_create2 (HD->dc, m, t_fft); /* complex */

  local Vec c_fft[m][m];
  bgy3d_vec_create2 (HD->dc, m, c_fft); /* complex */

  /* Solvent-solvent interaction  is a  sum of two  terms, short-range
     and long-range: */
  local Vec v_short[m][m];      /* real */
  bgy3d_vec_create2 (HD->da, m, v_short);

  local Vec v_long_fft[m][m];   /* complex */
  bgy3d_vec_create2 (HD->dc, m, v_long_fft);

  /* Get solvent-solvent site-site interactions: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      bgy3d_pair_potential (HD, solvent[i], solvent[j],
                            v_short[i][j], v_long_fft[i][j]);

  /* Prepare intra-molecular correlations: */
  local Vec w_fft[m][m];        /* diagonal will by NULL */
  bgy3d_omega_fft_create (HD, m, solvent, w_fft); /* creates them */

  /*
    For primary variable x there are  two ways to access the data: via
    the long Vec X and m * (m  + 1) / 2 shorter Vec x[m][m] aliased to
    the subsections of the longer one.
  */
  local Vec X = bgy3d_vec_pack_create2 (HD->da, m); /* long Vec */

  local Vec x[m][m];
  bgy3d_vec_aliases_create2 (X, m, x); /* aliases to subsections */

  /* Zero as intial guess  for x == t is not the  same as zero initial
     guess for x == c: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecSet (x[i][j], 0.0);

  /*
    Find an X such that dX  as returned by iterate_c2/t2 (&ctx, X, dX)
    is zero.  Cast is  there to  silence the mismatch  in the  type of
    first pointer argument: struct Ctx2* vs. void*:
  */
  {
    Ctx2 ctx =
      {
        .HD = HD,
        .m = m,
        .v_short = (void*) v_short,       /* in */
        .v_long_fft = (void*) v_long_fft, /* in */
        .y = (void*) y,                   /* work, t(c) or c(t) */
        .t_fft = (void*) t_fft,           /* work */
        .c_fft = (void*) c_fft,           /* work */
        .w_fft = (void*) w_fft,           /* in */
      };

    if (yes)
      bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_c2, X);
    else
      bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_t2, X);
  }

  /*
    This  saves   c??.bin  to  disk  to  be   used  in  solute/solvent
    calculations. The  other one  is not futhre  used, saved  just for
    info.
  */
  if (yes)
    {
      /* In  this branch  the primary  variable x  == c,  the indirect
         correlation t(c) is a dependent quantity y(x)! */
      bgy3d_vec_save2 ("c%d%d.bin", m, x);
      bgy3d_vec_save2 ("t%d%d.bin", m, y);
    }
  else
    {
      /* In  this branch  the  primary  variable x  ==  t, the  direct
         correlation c(t) is a dependent quantity y(x)! */
      bgy3d_vec_save2 ("t%d%d.bin", m, x);
      bgy3d_vec_save2 ("c%d%d.bin", m, y);
    }

  /*
    g  =  γ +  c  +  1,  store in  Vec  y.   The expression  is  valid
    irrespective of whether  (x, y) == (c,  γ) or (x, y) ==  (γ, c) as
    both appear as a sum x + y == γ + c here:
  */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        VecAXPY (y[i][j], 1.0, x[i][j]);
        VecShift (y[i][j], 1.0);
      }

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&t); */
  bgy3d_vec_destroy2 (m, t_fft);
  bgy3d_vec_destroy2 (m, c_fft);
  bgy3d_vec_destroy2 (m, v_short);
  bgy3d_vec_destroy2 (m, v_long_fft);

  /* Diagonal is NULL: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)
      bgy3d_vec_destroy (&w_fft[i][j]);

  bgy3d_vec_aliases_destroy2 (X, m, x);
  bgy3d_vec_pack_destroy2 (&X);

  bgy3d_state_destroy (HD);

  /* Dont forget to bgy3d_vec_destroy2() it! */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      g[i][j] = g[j][i] = y[i][j];

  bgy3d_vec_save2 ("g%d%d.bin", m, g);
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
  In k-space compute

   t = ρ c * h

  Here ρ  is a scalar  overall factor proportional to  solvent density
  and  should  include  convolution/FFT  renormalization in  order  to
  interprete the result as an FFT of t literally:
*/
static void compute_t1 (int m, real rho,
                        Vec c_fft[m][m], Vec h_fft[m], /* in */
                        Vec t_fft[m])                  /* out */
{
  /*
    fft(c)  *  fft(h).   Here   c  is  the  constant  (radial)  direct
    correlation  c2  of  the  pure solvent.   The  "convolution  star"
    corresponds to a matrix multiplication in the k-space.
  */
  complex pure fma (complex ti, complex cij, complex hj)
  {
    /* FMA stays for "fused multiply-add" */
    return ti + cij * hj;
  }

  /* For each solvent site ... */
  for (int i = 0; i < m; i++)
    {
      /* ... sum over solvent sites: */
      VecSet (t_fft[i], 0.0);
      for (int j = 0; j < m; j++)
        vec_fft_map3 (t_fft[i], /* argument aliasing! */
                      fma,
                      t_fft[i], c_fft[i][j], h_fft[j]);

      VecScale (t_fft[i], rho);
    }
}

/*
  HNC3d   iteration   for  a   fixed   direct   correlation  of   pure
  solvent. There are two cases:

  a)  Iteration  with   an  (almost  discontinous)  total  correlation
     function h as a primary variable

     h    ->  dh = h    - h
      in            out    in

  b) Iteration with the indirect correlation t as a primary variable

     t    ->  dt = t    - t
      in           out    in
*/
typedef struct Ctx1
{
  State *HD;
  int m;
  Vec *v;                       /* [m], real, fixed */
  Vec *y;                       /* [m], real, either h(t) or t(h) */
  Vec *h_fft, *t_fft;           /* [m], complex */
  Vec *c_fft;                   /* [m][m], complex, fixed */
} Ctx1;


static void iterate_h1 (Ctx1 *ctx, Vec H, Vec dH)
{
  /* alias of the right shape: */
  const int m = ctx->m;
  Vec (*c_fft)[m] = (void*) ctx->c_fft;
  Vec *t = (void*) ctx->y;      /* [m], here y = t(h) */

  const ProblemData *PD = ctx->HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real N3 = PD->N[0] * PD->N[1] * PD->N[2];

  /* Establish aliases to the subsections of the long Vec H and dH: */
  local Vec h[m], dh[m];
  bgy3d_vec_aliases_create1 (H, m, h);
  bgy3d_vec_aliases_create1 (dH, m, dh);

  /* fft(h).  Here h is the 3d  unknown hole density h1 of the solvent
     sites. */
  for (int i = 0; i < m; i++)
    MatMult (ctx->HD->fft_mat, h[i], ctx->h_fft[i]);

  /*
    fft(c)  *  fft(h).   Here   c  is  the  constant  (radial)  direct
    correlation  c2  of  the  pure solvent.   The  "convolution  star"
    corresponds to a matrix multiplication in the k-space.

    For  each solvent  site sum  over solvent  sites. Let  the overall
    scale include  a factor for  inverse FFT and  convolution integral
    right away:
  */
  compute_t1 (m, rho * h3 / N3, c_fft, ctx->h_fft, ctx->t_fft);

  /* t = fft^-1 (fft(c) * fft(h)). Here t is 3d t1. */
  for (int i = 0; i < m; i++)
    MatMultTranspose (ctx->HD->fft_mat, ctx->t_fft[i], t[i]);

  /*
    The new candidate for the total correlation

      h ~ exp (-βv + t) - 1

    See comments to compute_h() for the meaning of "~".
  */
  for (int i = 0; i < m; i++)
    compute_h (beta, ctx->v[i], t[i], dh[i]);

  /*
    This  destroys the  aliases,  but  does not  free  the memory,  of
    course. The actuall data is owned by Vec H and Vec dH. From now on
    one may access H and dH directly again.
  */
  bgy3d_vec_aliases_destroy1 (H, m, h);
  bgy3d_vec_aliases_destroy1 (dH, m, dh);

  /*
    dh := h    - h
           out    in
  */
  VecAXPY (dH, -1.0, H);
}


static void iterate_t1 (Ctx1 *ctx, Vec T, Vec dT)
{
  /* alias of the right shape: */
  const int m = ctx->m;
  Vec (*c_fft)[m] = (void*) ctx->c_fft;
  Vec *h = (void*) ctx->y;      /* [m], here y = h(t) */

  const ProblemData *PD = ctx->HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real N3 = PD->N[0] * PD->N[1] * PD->N[2];

  /* Establish aliases to the subsections of the long Vec T and dT: */
  local Vec t[m], dt[m];
  bgy3d_vec_aliases_create1 (T, m, t);
  bgy3d_vec_aliases_create1 (dT, m, dt);

  /*
    The new candidate for the total correlation

      h ~ exp (-βv + t) - 1

    See comments to compute_h() for the meaning of "~".
  */
  for (int i = 0; i < m; i++)
    {
      compute_h (beta, ctx->v[i], t[i], h[i]);

      /* fft(h).   Here h is  the 3d  unknown hole  density h1  of the
         solvent sites. */
      MatMult (ctx->HD->fft_mat, h[i], ctx->h_fft[i]);
    }

  /*
    fft(c)  *  fft(h).   Here   c  is  the  constant  (radial)  direct
    correlation  c2  of  the  pure solvent.   The  "convolution  star"
    corresponds to a matrix multiplication in the k-space.

    For  each solvent  site sum  over solvent  sites. Let  the overall
    scale include  a factor for  inverse FFT and  convolution integral
    right away:
  */
  compute_t1 (m, rho * h3 / N3, c_fft, ctx->h_fft, ctx->t_fft);

  /* t = fft^-1 (fft(c) * fft(h)). Here t is 3d t1. */
  for (int i = 0; i < m; i++)
    MatMultTranspose (ctx->HD->fft_mat, ctx->t_fft[i], dt[i]);

  /*
    This  destroys the  aliases,  but  does not  free  the memory,  of
    course. The actuall data is owned by Vec T and Vec dT. From now on
    one may access T and dT directly again.
  */
  bgy3d_vec_aliases_destroy1 (T, m, t);
  bgy3d_vec_aliases_destroy1 (dT, m, dt);

  /*
    dt := t    - t
           out    in
  */
  VecAXPY (dT, -1.0, T);
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


/*
  Solving  HNC3d  equations.  Direct  correlation  c  of pure  solvent
  appears as a fixed input here. A  primary variable is either h = g -
  1 or t related to h by one of the closure relations.

  At least for  water in OW solvent (water-like  model, just a neutral
  oxygen with LJ parameters of water) Newton iteration over t ~ log [g
  exp(βv)] converges, whereas iteration over h  = g - 1 does not.  The
  problem with h  as a primary variable cannot be  severe as the Jager
  solver  manages  that. On  the  other  hand h  =  g  -  1 is  nearly
  discontinuous at the cavity border and formally constrained h >= -1.
  Intuitively  this does  not make  the  task of  the fxipoint  solver
  easier.
*/
void hnc3d_solute_solve (const ProblemData *PD,
                         const int m, const Site solvent[m],
                         const int n, const Site solute[n],
                         Vec g[m])
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  /* Flip this to switch between h and t as primary variables: true =>
     h, false => t: */
  const bool yes = false;

  if (yes)
    PetscPrintf (PETSC_COMM_WORLD, "(iterations for h)\n");
  else
    PetscPrintf (PETSC_COMM_WORLD, "(iterations for γ)\n");

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */

  /*
    This will be a functional y(x) of primary variable, either y(x) ==
    t(h)  or  vice  versa, y(x)  =  h(t).  Depending  on the  type  of
    iteration.
  */
  Vec y[m];                     /* FIXME: not deallocated */
  bgy3d_vec_create1 (HD->da, m, y);

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

  /*
    For primary variable x there are  two ways to access the data: via
    the long Vec  X and m shorter Vec x[m]  aliased to the subsections
    of the longer one.
  */
  local Vec X = bgy3d_vec_pack_create1 (HD->da, m); /* long Vec */

  local Vec x[m];
  bgy3d_vec_aliases_create1 (X, m, x); /* aliases to subsections */

  /* Zero as intial guess for x == t is (almost?) the same as exp(-βv)
     - 1 initial guess for x == h: */
  if (yes)
    for (int i = 0; i < m; i++)
      compute_h0 (HD->PD->beta, v[i], x[i]);
  else
    for (int i = 0; i < m; i++)
      VecSet (x[i], 0.0);

  /*
    Find a  t such that  dt as returned  by iterate_t1 (HD, t,  dt) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: Ctx1* vs. void*:
  */
  {
    /* Work area for iterate_h1(): */
    Ctx1 ctx =
      {
        .HD = HD,
        .m = m,
        .v = (void*) v,
        .y = (void*) y,         /* h(t) or t(h) */
        .c_fft = (void*) c_fft, /* pair quantitity */
        .h_fft = (void*) h_fft,
        .t_fft = (void*) t_fft,
      };

    if (yes)
      bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_h1, X);
    else
      bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_t1, X);
  }

  /* Vec  y[]  will  be  returned  to  the caller,  Vec  x[]  will  be
     destroyed: */
  if (yes)
    for (int i = 0; i < m; i++)
      VecCopy (x[i], y[i]);
  /*
    else nothing, but read the note ...

    FIXME: by ignoring the value of Vec x == t and using Vec y == h(t)
    instead, we  assume that the  two are consistent  as in y  = y(x).
    Hopefully the SNES solver does not return "untested" solutions.  A
    call to iterate_t1/h1() to "test"  a solution would update Vec y =
    y(x), namely.
  */

  /* g := h + 1 */
  for (int i = 0; i < m; i++)
    VecShift (y[i], 1.0);

  /* free stuff */
  /* Delegated to the caller: bgy3d_vec_destroy (&h) */
  bgy3d_vec_destroy1 (m, h_fft);
  bgy3d_vec_destroy1 (m, t_fft);
  bgy3d_vec_destroy1 (m, v);

  bgy3d_vec_aliases_destroy1 (X, m, x);
  bgy3d_vec_pack_destroy1 (&X);

  /* This should be the only pair quantity: */
  bgy3d_vec_destroy2 (m, c_fft);

  bgy3d_state_destroy (HD);

  /* Caller is supposed to destroy it! */
  for (int i = 0; i < m; i++)
    g[i] = y[i];

  bgy3d_vec_save1 ("g%d.bin", m, g);
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
