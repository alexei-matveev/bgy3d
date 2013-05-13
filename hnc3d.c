/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: hnc3d.c,v 1.13 2006-12-14 17:35:38 jager Exp $ */
/*==========================================================*/

/*
  For simple liquids the matrix OZ equation reads

    h = c * [1 + ρ h]

  where the expression in the  square brackets is susceptibility χ = 1
  +  ρh of  the mixture.   The matrix  ρ is  diagonal.  Only  when the
  density of all sites is the same the equation simplifies to the more
  recognizable form

    h = c + ρ c * h.

  However, consider a solute/solvent  mixture in the limit of infinite
  dilution where the  density of the solute is  vanishingly small. The
  corresponding blocks of the matrix equation read:

    h   = c   [1 + ρ  h  ] + c   ρ  h                              (1)
     vv    vv       v  vv     vu  u  uv

  The second  term in  Eq. (1)  vanishes and the  OZ equation  for the
  solvent   structure  becomes   independent  of   the   solute.   The
  solute-solute  correlations are,  of course,  affected by  the dense
  medium even in the infinite dilution limit:

    h   = c   [1 + ρ  h  ] + c   ρ  h                              (2)
     uu    uu       u  uu     uv  v  vu

  Now the remaining two equations

    h   = c   [1 + ρ  h  ] + c   ρ  h                              (3)
     uv    uv       v  vv     uu  u  uv

  and

    h   = c   [1 + ρ  h  ] + c   ρ  h                              (4)
     vu    vu       u  uu     vv  v  vu

  though necessarily  redundant seem to differ. By  removing the terms
  proportional  to  the solute  density  one  gets  two different  but
  equivalent relations:

    h   = c   + c   ρ  h                                           (5)
     uv    uv    uv  v  vv

  and

    h   = c   + c   ρ  h                                           (6)
     vu    vu    vv  v  vu

  For  practical application  these two  equations  may be  used as  a
  relation between the  solute-solvent indirect correlation t =  h - c
  and a another solute-solvent property as a complement to the closure
  relation.  Both equation use the solvent structure as a parameter:

    t   = c   ρ  h                                                 (7)
     uv    uv  v  vv

  or

    t   = c   ρ  h                                                 (8)
     vu    vv  v  vu

  In  the first  case the  structure of  the pure  solvent  (or rather
  infinitely  diluted mixture) is  represented by  the solvent-solvent
  total  correlation  h,  whereas  in  the  second  case  the  solvent
  structure is  represented by the  solvent-solvent direct correlation
  c. The second  of these two forms, Eq.  (8), was originally proposed
  for  use  in 3d  HNC  equations  [1].  However other  (more  recent)
  literature  appears  to prefer  to  start  with  a formulation  that
  resembles more  Eq.  (7) e.g.  Ref. [2].   This alternative approach
  assumes that the structure of the pure solvent is represented by the
  solvent  susceptibility χ  that  mediates the  relation between  the
  total and direct solute-solvent correlations:

    h   = c   χ
     uv    uv  vv

  with

    χ   = 1 + ρ  h
     uv        v  vv

  which is equivalent  to Eq. (7).

  The discussion  above deals  with a special  case of a  more general
  molecular RISM equation

    h   =  ω  * c   * [ω + ρ  h  ]
     uv     u    uv     v   v  vv

  Where  again the  term is  square brackets  is the  property  of the
  molecular solvent fixed by  the assumption of the infinite dilution,
  and is a generalization  of the solvent susceptibility for molecular
  solvents:

    χ   =  ω  + ρ  h
     vv     v    v  vv

  [1] Numerical solution of the hypernetted chain equation for a
      solute of arbitrary geometry in three dimensions, Dmitrii Beglov
      and Benoît Roux , J.  Chem.  Phys.  103, 360 (1995);
      http://dx.doi.org/10.1063/1.469602

  [2] Three-Dimensional Molecular Theory of Solvation Coupled with
      Molecular Dynamics in Amber, Tyler Luchko, Sergey Gusarov,
      Daniel R. Roe, Carlos Simmerling, David A. Case , Jack Tuszynski
      and Andriy Kovalenko, J. Chem. Theory Comput., 2010, 6 (3), pp
      607–624 http://dx.doi.org/10.1021/ct900460m
*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-fftw.h"         /* bgy3d_fft_mat_create() */
#include "bgy3d-vec.h"          /* vec_create() */
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

  /* Here ij  and ji are aliased.   We ask to put  real* into complex*
     slots: */
  vec_get_array2 (m, c_fft, (void*) c_fft_);
  vec_get_array2 (m, t_fft, (void*) t_fft_);

  /* Diagonals are NULL, vec_get_array() passes them through: */
  vec_get_array2 (m, w_fft, (void*) w_fft_);

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

  /* Here vec_restore_array2() expects array  of real*, we offer array
     of complex* instead: */
  vec_restore_array2 (m, c_fft, (void*) c_fft_);
  vec_restore_array2 (m, t_fft, (void*) t_fft_);

  /* The case with NULL on the diagonal is handled gracefully: */
  vec_restore_array2 (m, w_fft, (void*) w_fft_);
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
  There  were  historically  two  types  of HNC  iterations  for  pure
  solvent:

  a) HNC iteration for an indirect correlation t = h - c

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
     of c and is stored in y = t(c).

  The case (a)  with indirect correlation t as  a primary variable won
  over time. Maintaining two  cases appeared infeasible at some point.
*/
typedef struct Ctx2
{
  State *HD;
  int m;
  Vec *v_short;              /* [m][m], real, center, const */
  Vec *v_long_fft;           /* [m][m], complex, center, const */
  Vec *y;                    /* [m][m], real, y = c(t), not t(c) */
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
  vec_aliases_create2 (T, m, t);
  vec_aliases_create2 (dT, m, dt);

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
  vec_aliases_destroy2 (T, m, t);
  vec_aliases_destroy2 (dT, m, dt);

  VecAXPY (dT, -1.0, T);
}

/*
 * h(r) = c(r) + γ(r)
 * μ := ½h²(r) - c(r) - ½h(r)c(r)
 */
static void compute_mu (Vec c, Vec t, Vec mu)
{
  real pure f (real c, real t)
  {
    /* h = c + γ */
    return 0.5 * (c + t) * (c + t) - c - 0.5 * (c + t) * c;
  }
  vec_map2 (mu, f, c, t);
}

/*
  Returns the β-scaled "density" of the chemical potential, βμ(r).
  To  get the  excess  chemical potential  integrate  it over  the
  volume and divide by β, cf.:

  βμ = 4πρ ∫ [½h²(r) - c(r) - ½h(r)c(r)] r²dr

  Here we pass c(r) and γ(r) = h(r) - c(r)

  Volume integral in cartesian grid is actually:

  Vol(D) = ∫∫∫dxdydz
            D
*/
static void chempot_density (int m,
                            Vec c[m][m], Vec t[m][m], /* in */
                            Vec mu)                   /* out */
{
  /* increment for all solvent sites */
  local Vec dmu = vec_duplicate (mu);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
    {
      compute_mu (c[i][j], t[i][j], dmu);
      VecAXPY (mu, 1.0, dmu);
    }

  vec_destroy (&dmu);
}

/*
  Solving for indirect correlation t = h - c and thus, also for direct
  correlation c  and other quantities  of HNC equation.   The indirect
  correlation  t appears  as  a  primary variable  x  here. All  other
  unknowns, including the direction  correlation c, are functionals of
  that.

  Historically  this  code supported  another  mode  where the  direct
  correlation  c served as  a primary  variable x.   In that  case the
  indirect correlation t was a functional of that.
*/
void hnc3d_solvent_solve (const ProblemData *PD,
                          int m, const Site solvent[m],
                          Vec g[m][m])
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */

  PetscPrintf (PETSC_COMM_WORLD, "(iterations for γ)\n");

  /*
    This  will be  a functional  y(x) of  primary variable  x,  y(x) =
    c(t). Earlier versions supported direct correlation c as a primary
    variable  so  that  in that  case  one  had  y(x)  == t(c)  as  an
    alternative.
  */
  Vec y[m][m];                  /* FIXME: not deallocated */
  vec_create2 (HD->da, m, y);

  local Vec t_fft[m][m];
  vec_create2 (HD->dc, m, t_fft); /* complex */

  local Vec c_fft[m][m];
  vec_create2 (HD->dc, m, c_fft); /* complex */

  /* Solvent-solvent interaction  is a  sum of two  terms, short-range
     and long-range: */
  local Vec v_short[m][m];      /* real */
  vec_create2 (HD->da, m, v_short);

  local Vec v_long_fft[m][m];   /* complex */
  vec_create2 (HD->dc, m, v_long_fft);

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
  local Vec X = vec_pack_create2 (HD->da, m); /* long Vec */

  local Vec x[m][m];
  vec_aliases_create2 (X, m, x); /* aliases to subsections */

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
        .y = (void*) y,                   /* work, c(t), not t(c) */
        .t_fft = (void*) t_fft,           /* work */
        .c_fft = (void*) c_fft,           /* work */
        .w_fft = (void*) w_fft,           /* in */
      };

    bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_t2, X);
  }

  /*
    This saves c??-fft.bin (including  the long-range part) to disk to
    be  used in  solute/solvent calculations.   The approach  with the
    primary variable  x ==  t won over  time.  The  direct correlation
    c(t) is a dependent quantity y(x)!
  */
  {
    Vec (*c)[m] = y;            /* read-only alias */
    const real dV = PD->h[0] * PD->h[1] * PD->h[2];

    /*
      The real-space rep of the short-range direct correlation is only
      usefull for debug and visualization:

      bgy3d_vec_save2 ("c%d%d.bin", m, c);
    */

    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i; j++)
        {
          MatMult (HD->fft_mat, c[i][j], c_fft[i][j]);
          VecScale (c_fft[i][j], dV);

          /*
            The real-space representation encodes only the short-range
            part  of the  direct  corrlation. The  (fixed) long  range
            contribution is added here:

              C := C  - βV
                    S     L
          */
          VecAXPY (c_fft[i][j], -PD->beta, v_long_fft[i][j]);

          /* Translate the  distribution to the grid  corner.  This is
             what one expects in convolution integrals: */
          bgy3d_vec_fft_trans (HD->dc, HD->PD->N, c_fft[i][j]);
        }

    /* The Fourier  rep is what is actually  read in solvent_kernel(),
       see below: */
    bgy3d_vec_save2 ("c%d%d-fft.bin", m, c_fft);
  }

  /* chemical potential */
  {
    const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
    local Vec mu_dens = vec_create (HD->da);

    /* x == t; y == c */
    chempot_density (m, y, x, mu_dens);

    /* Volume integral scaled by a factor: */
    const real mu = PD->rho * vec_sum (mu_dens) * h3 / PD->beta;

    PetscPrintf (PETSC_COMM_WORLD, " mu = %f\n", mu);

    bgy3d_vec_save ("mu_dens.bin", mu_dens);
    vec_destroy (&mu_dens);
  }


  /*
    h = c  + t, store in Vec y.  The  expression is valid irrespective
    of whether (x, y) == (c, t) or  (x, y) == (t, c) as both appear as
    a sum x + y == t + c here:
  */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecAXPY (y[i][j], 1.0, x[i][j]);

  /* g = 1 + h, store in Vec y for output: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecShift (y[i][j], 1.0);

  /* free stuff */
  /* Delegated to the caller: vec_destroy (&t); */
  vec_destroy2 (m, t_fft);
  vec_destroy2 (m, c_fft);
  vec_destroy2 (m, v_short);
  vec_destroy2 (m, v_long_fft);

  /* Diagonal is NULL: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)
      vec_destroy (&w_fft[i][j]);

  vec_aliases_destroy2 (X, m, x);
  vec_pack_destroy2 (&X);

  bgy3d_state_destroy (HD);

  /* Dont forget to vec_destroy2() it! */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      g[i][j] = g[j][i] = y[i][j];

  bgy3d_vec_save2 ("g%d%d.bin", m, g);
}


/*
  Solving for  g of HNC  equation with non-linear solver  specified by
  the --snes-solver flag.
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
  return g[0][0];               /* vec_destroy (&) it! */
}


/*
  In k-space compute

   t   = ρ c   *  h
    vu      vv     vu

  Compare  to Eq.  (8) in  the  header comments.  Here ρ  is a  scalar
  overall factor  equal to solvent density if  the convolution theorem
  is factorless.  The caller may  (ab)use this factor to further scale
  the result.

  Note that the solvent-solvent  site-site correlation c is eventually
  long-range.  In Fourier rep that would  mean c(k) is singular at k =
  0. Here  there are no  special precautions to treat  the assymptotic
  interactions of the (so far neutral) solute species with the charged
  solvent sites. Note that in 3d solute/solvent case we do not operate
  with  site-site   distributions/potentials,  but  rather   with  the
  "solvent site"-"compound solute" quantities. If the physics works as
  expected the assymptotic decay  of such potentials and distributions
  is "fast"  for neutral solutes. FIXME:  this is not  true anymore if
  the solute is charged.
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
  solvent. There were two cases at some point:

  a) Iteration with the indirect correlation t as a primary variable

     t    ->  dt = t    - t
      in            out    in

  b)  Iteration  with   an  (almost  discontinous)  total  correlation
     function h as a primary variable

     h    ->  dh = h    - h
      in            out    in

  At some point  it appeared impractical to maintain  two branches, so
  that the case (a) with  indirect correlation t as a primary variable
  won.
*/
typedef struct Ctx1
{
  State *HD;
  int m;
  Vec *v;                       /* [m], real, fixed */
  Vec *y;                       /* [m], real, y = h(t), not t(h) */
  Vec *h_fft, *t_fft;           /* [m], complex */
  Vec *c_fft;                   /* [m][m], complex, fixed */
} Ctx1;


static void iterate_t1 (Ctx1 *ctx, Vec T, Vec dT)
{
  /* alias of the right shape: */
  const int m = ctx->m;
  Vec (*c_fft)[m] = (void*) ctx->c_fft;
  Vec *h = (void*) ctx->y;      /* [m], here y = h(t) */

  const ProblemData *PD = ctx->HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real N3 = PD->N[0] * PD->N[1] * PD->N[2];

  /* Establish aliases to the subsections of the long Vec T and dT: */
  local Vec t[m], dt[m];
  vec_aliases_create1 (T, m, t);
  vec_aliases_create1 (dT, m, dt);

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
    scale include a factor for inverse FFT right away:
  */
  compute_t1 (m, rho / N3, c_fft, ctx->h_fft, ctx->t_fft);

  /* t = fft^-1 (fft(c) * fft(h)). Here t is 3d t1. */
  for (int i = 0; i < m; i++)
    MatMultTranspose (ctx->HD->fft_mat, ctx->t_fft[i], dt[i]);

  /*
    This  destroys the  aliases,  but  does not  free  the memory,  of
    course. The actuall data is owned by Vec T and Vec dT. From now on
    one may access T and dT directly again.
  */
  vec_aliases_destroy1 (T, m, t);
  vec_aliases_destroy1 (dT, m, dt);

  /*
    dt := t    - t
           out    in
  */
  VecAXPY (dT, -1.0, T);
}


/* Reads  c_fft[][]  as previousely  written  by  solvent solver.  See
   hnc3d_solvent_solve() above. */
static void solvent_kernel (State *HD, int m, Vec c_fft[m][m])
{
  if (bgy3d_getopt_test ("--from-radial-g2")) /* FIXME: better name? */
    {
      /*
        Load  radial  direct   correlation  from  text  file.   FIXME:
        representing long-range on  a real-space grid is intrinsically
        broken!
      */
      const real dV = HD->PD->h[0] * HD->PD->h[1] * HD->PD->h[2];

      local Vec c[m][m];
      vec_create2 (HD->da, m, c);

      bgy3d_vec_read_radial2 (HD->da, HD->PD, "c%d%d.txt", m, c);

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          {
            MatMult (HD->fft_mat, c[i][j], c_fft[i][j]);
            VecScale (c_fft[i][j], dV);

            /* Translate the distribution to  the grid corner. This is
               what one expects in convolution integrals. */
            bgy3d_vec_fft_trans (HD->dc, HD->PD->N, c_fft[i][j]);
          }
      vec_destroy2 (m, c);
    }
  else
    bgy3d_vec_read2 ("c%d%d-fft.bin", m, c_fft); /* ready for use as is */

}


/*
  Solving  HNC3d  equations.  Direct  correlation  c  of pure  solvent
  appears  as a  fixed  input  here. A  primary  variable is  indirect
  correlation t related to h by one of the closure relations.

  Historically  another   branch  of   the  code  treated   the  total
  correlation  h as  a primary  variable.  At  least for  water  in OW
  solvent (water-like model, just  a neutral oxygen with LJ parameters
  of  water) Newton  iteration over  t  ~ log  [g exp(βv)]  converges,
  whereas iteration over h = g -  1 does not.  The problem with h as a
  primary  variable  cannot be  severe  as  the  Jager solver  manages
  that. On  the other hand h  = g -  1 is nearly discontinuous  at the
  cavity border and formally constrained by h >= -1.  Intuitively this
  does  not make  the  task of  the  fixpoint solver  easier. So  that
  alternative has been abandoned in favor of indirect correlation t as
  a primary variable.
*/
void hnc3d_solute_solve (const ProblemData *PD,
                         const int m, const Site solvent[m],
                         const int n, const Site solute[n],
                         Vec g[m])
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  /*
    This  will be  a functional  y(x) of  primary variable  x,  y(x) =
    c(t). Earlier versions supported direct correlation c as a primary
    variable  so  that  in that  case  one  had  y(x)  == t(c)  as  an
    alternative.
  */
  PetscPrintf (PETSC_COMM_WORLD, "(iterations for γ)\n");

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */

  /*
    This will be a functional y(x) of primary variable, either y(x) ==
    t(h)  or  vice  versa, y(x)  =  h(t).  Depending  on the  type  of
    iteration.
  */
  Vec y[m];                     /* FIXME: not deallocated */
  vec_create1 (HD->da, m, y);

  local Vec h_fft[m];
  vec_create1 (HD->dc, m, h_fft); /* complex */

  local Vec t_fft[m];
  vec_create1 (HD->dc, m, t_fft); /* complex */

  local Vec v[m];
  vec_create1 (HD->da, m, v); /* solute-solvent interaction */

  /* This should be the only pair quantity: */
  local Vec c_fft[m][m];
  vec_create2 (HD->dc, m, c_fft);        /* complex */

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
  local Vec X = vec_pack_create1 (HD->da, m); /* long Vec */

  local Vec x[m];
  vec_aliases_create1 (X, m, x); /* aliases to subsections */

  /* Zero as intial guess for x == t is (almost?) the same as exp(-βv)
     - 1 initial guess for x == h: */
  for (int i = 0; i < m; i++)
    VecSet (x[i], 0.0);

  /*
    Find a  t such that  dt as returned  by iterate_t1 (HD, t,  dt) is
    zero. Cast is  there to silence the mismatch in  the type of first
    pointer argument: Ctx1* vs. void*:
  */
  {
    /* Work area for iterate_t1(): */
    Ctx1 ctx =
      {
        .HD = HD,
        .m = m,
        .v = (void*) v,
        .y = (void*) y,         /* h(t), not t(h) */
        .c_fft = (void*) c_fft, /* pair quantitity */
        .h_fft = (void*) h_fft,
        .t_fft = (void*) t_fft,
      };

    bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_t1, X);
  }

  /*
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
  /* Delegated to the caller: vec_destroy (&h) */
  vec_destroy1 (m, h_fft);
  vec_destroy1 (m, t_fft);
  vec_destroy1 (m, v);

  vec_aliases_destroy1 (X, m, x);
  vec_pack_destroy1 (&X);

  /* This should be the only pair quantity: */
  vec_destroy2 (m, c_fft);

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
