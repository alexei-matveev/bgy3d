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

    h   = c   χ                                                    (9)
     uv    uv  vv

  with

    χ   = 1 + ρ  h                                                 (10)
     uv        v  vv

  which is  equivalent to Eq. (7). Note  that χ, as it  stands, is not
  symmetric, unless ρ is the same for all sites. FIXME: the code makes
  use of this assumption!

  The discussion  above deals  with a special  case of a  more general
  molecular RISM equation

    h   =  ω  * c   * [ω + ρ  h  ]                                (11)
     uv     u    uv     v   v  vv

  Where  again the  term is  square brackets  is the  property  of the
  molecular solvent fixed by  the assumption of the infinite dilution,
  and is a generalization  of the solvent susceptibility for molecular
  solvents:

    χ   =  ω  + ρ  h                                              (12)
     vv     v    v  vv

  With this definition, the expression for the indirect solute-solvent
  correlation

    t   = h   - c   = c   (χ   - 1)                               (13)
     uv    uv    uv    uv   vv

  does  not simplify  to  Eq. (7)  directly  though still  has a  very
  similar structure. Again,  only if ρ is a scalar  one can assume the
  solvent susceptibility χ is symmetric in site indices.

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
#include "rism.h"               /* rism_solvent() */
#include "bgy3d-potential.h"    /* info() */
#include "bgy3d-solvents.h"      /* bgy3d_sites_show() */
#include <math.h>               /* expm1() */
#include "hnc3d.h"


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


static void compute_c (ClosureEnum closure, real beta, Vec v, Vec t, Vec c)
{
  switch (closure)
    {
    case CLOSURE_HNC:
      compute_c_HNC (beta, v, t, c);
      break;
    case CLOSURE_KH:
      compute_c_KH (beta, v, t, c);
      break;
    case CLOSURE_PY:
      compute_c_PY (beta, v, t, c);
      break;
    }
  /* No  default,  let the  compiler  warn if  we  do  not handle  all
     cases. */
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
      in            out    in

     where other  intermediates are considered a function  of that. In
     particular  in this case  another important  intermediate is  y =
     c(t).

  b) HNC iteration for a direct correlation where the total correlation
     is considered a function of that:

     c    ->  dc = c    - c
      in            out    in

     In this case the indirect  correlation γ is considered a function
     of c and is stored in y = t(c).

  The case (a)  with indirect correlation t as  a primary variable won
  over time. Maintaining two  cases appeared infeasible at some point.
*/
typedef struct Ctx2
{
  State *HD;
  int m;
  Vec *v_short;       /* [m][m], in, real, center */
  Vec *v_long_fft;    /* [m][m], in, complex, center */
  Vec *y;             /* [m][m], work, real, y = c(t), not t(c) */
  Vec *t_fft, *c_fft; /* [m][m], work, complex */
  Vec *w_fft;         /* [m][m], in, complex, corner, NULL diagonal */
} Ctx2;


static void iterate_t2 (Ctx2 *ctx, Vec T, Vec dT)
{
  const int m = ctx->m;

  /* aliases of the correct [m][m] shape: */
  Vec (*c_fft)[m] = (void*) ctx->c_fft;
  Vec (*t_fft)[m] = (void*) ctx->t_fft;
  Vec (*w_fft)[m] = (void*) ctx->w_fft; /* corner */
  Vec (*c)[m] = (void*) ctx->y;         /* y = c(t) here */
  Vec (*v_short)[m] = (void*) ctx->v_short;
  Vec (*v_long_fft)[m] = (void*) ctx->v_long_fft;

  const State *HD = ctx->HD;
  const ProblemData *PD = HD->PD;
  const real rho = PD->rho;
  const real beta = PD->beta;
  const real L3 = PD->L[0] * PD->L[1] * PD->L[2];
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  /* Establish aliases to the subsections of the long Vec T and dT: */
  local Vec t[m][m], dt[m][m];
  vec_aliases_create2 (T, m, t);
  vec_aliases_create2 (dT, m, dt);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        compute_c (PD->closure, beta, v_short[i][j], t[i][j], c[i][j]);

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

        VecScale (dt[i][j], 1.0/L3);
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
         2
    μ = h / 2 - c  - h * (c + c ) / 2
                 S         S   L

  This is  a subexpression in  the summation over all  site-site pairs
  where it is assumed that the long-range contribution in

    - Σ   c
       uv  uv

  cancels out. This is the case when

    c      = -βQ q v
     L,uv       u v L

  and either  solute charges Q or  the solvent charges q  sum to zero.
  Naturally, with  an exact  arithmetics one does  not need  to handle
  short- and long-range correlations separately in this case.
*/
static void compute_mu (real thresh, Vec h, Vec cs, Vec cl, Vec mu)
{
  /*
    The h² term  contributes conditionally. Eventually, only depletion
    regions (h < 0) contribute  (KH).  Threshold is supposed to be 0.0
    for KH functional (depletion regions contribute), anywhere between
    1  and  +∞  for GF  functional  (no  such  term)  and -∞  for  HNC
    functional (contributes unconditionally):
  */
  real pure f (real h, real cs, real cl)
  {
    const real h2 = (-h > thresh ? h * h : 0.0);
    return h2 / 2 - cs - h * (cs + cl) / 2;
  }
  vec_map3 (mu, f, h, cs, cl);
}


/*
  Returns the  β-scaled density of the chemical  potential, βμ(r).  To
  get the excess  chemical potential integrate it over  the volume and
  divide by β, cf.:

    βμ = ρ ∫ [½h²(r) - c(r) - ½h(r)c(r)] d³r

  Here  we  pass  h(r)  and  c(r)  and long-range  part  of  c(r)  for
  solute-solvent pair.  The solute  might be the  same as  solvent, of
  course.

  Volume integral in cartesian grid is actually:

    Vol(D) = ∫∫∫dxdydz
              D

  The shape of arrays passed here may be arbitrary, though in practice
  it is either m x m or 1 x m:
*/
static void chempot_density (real thresh, int n, int m,
                             Vec h[n][m],  /* in */
                             Vec cs[n][m], /* in */
                             Vec cl[n][m], /* in, long range */
                             Vec mu)       /* out */
{
  /* increment for all solvent sites */
  local Vec dmu = vec_duplicate (mu);

  /* Clear accumulator: */
  VecSet (mu, 0.0);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      {
        compute_mu (thresh, h[i][j], cs[i][j], cl[i][j], dmu);

        VecAXPY (mu, 1.0, dmu);
      }

  vec_destroy (&dmu);
}


/*
  Interface to get chemical potential of solvent-solvent pair from Vec
  h,  c, and cl,  return the  value of  chemical potential.

  The   argument   "closure"   choses   between   HNC,   KH   and   GF
  functionals. The State  struct also holds a setting  for the closure
  --- that one is  ignored to allow computing any  kind of functional.
*/
static real chempot (const State *HD, ClosureEnum closure, int n, int m,
                     Vec h[n][m], Vec c[n][m], Vec cl[n][m]) /* in */
{
  const ProblemData *PD = HD->PD;
  const real beta = PD->beta;
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  real thresh;
  switch (closure)              /* not PD->closure! */
    {
    case CLOSURE_KH:
      /* The  h²  term  contributes  only  in  the  depletion  regions
         (KH): */
      thresh = 0.0;
      break;
    case CLOSURE_HNC:
      /* The h² term contributes unconditionally (HNC): */
      thresh = - DBL_MAX;
      break;
    default:
      /* There is no h² term otherwise (GF): */
      thresh = DBL_MAX;
      break;
    }

  /* Vector for chemical potential density */
  local Vec mu_dens = vec_create (HD->da);

  /* Get β-scaled chemical potential density */
  chempot_density (thresh, n, m, h, c, cl, mu_dens);

  /* Volume integral scaled by a factor: */
  const real mu = PD->rho * vec_sum (mu_dens) * h3 / beta;

  vec_destroy (&mu_dens);

  return mu;
}


/* Prints chemical  potentials on  tty. The default  is marked  with a
   star. */
static void
print_chempot (const State *HD, int n, int m,
               Vec h[n][m], Vec c[n][m], Vec cl[n][m]) /* in */
{
  /* FIXME:  as  implemented, for  method  ==  PY  the default  energy
     functional is GF: */
  const ClosureEnum methods[3] = {CLOSURE_HNC, CLOSURE_KH, CLOSURE_PY};
  const char *names[3] = {"HNC", "KH", "GF"};
  real mu[3];

  /* Silent computing: */
  for (int i = 0; i < 3; i++)
    mu[i] = chempot (HD, methods[i], n, m, h, c, cl);

  /* Printing only: */
  PetscPrintf (PETSC_COMM_WORLD,
               " # Chemical potentials, default is marked with *:\n");
  for (int i = 0; i < 3; i++)
    PetscPrintf (PETSC_COMM_WORLD, " # mu = %f kcal (%s)%s\n", mu[i], names[i],
                 ((methods[i] == HD->PD->closure)? "*" : ""));
}


/*
  Solving for indirect correlation t = h - c and thus, also for direct
  correlation c  and other quantities  of HNC equation.   The indirect
  correlation  t  appears  as  a  primary  variable  here.  All  other
  unknowns, including the direction  correlation c, are functionals of
  that.

  Historically  this  code supported  another  mode  where the  direct
  correlation  c served  as  a  primary variable.   In  that case  the
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
    For primary variable t there  will be two exclusive ways to access
    the data:  via the  long Vec  T or  m * (m  + 1)  / 2  shorter Vec
    t[m][m] aliased to the subsections of the longer one.
  */
  local Vec T = vec_pack_create2 (HD->da, m); /* long Vec */

#ifndef WITH_FORTRAN
  /* Zero as intial guess for t  is not the same as zero initial guess
     for c: */
  VecSet (T, 0.0);
#else
  {
    /*
      Problem data for use in 1d-RISM code. If you do not upscale, the
      rmax = max (PD->L) / 2 will be too low for interpolation:
    */
    ProblemData pd = rism_upscale (PD);

    const int nrad = rism_nrad (&pd);
    const real rmax = rism_rmax (&pd);

    real t_rad[m][m][nrad];

    /* 1d-RISM calculation.  Dont need susceptibility, need t(r): */
    rism_solvent (&pd, m, solvent, t_rad, NULL);

    {
      local Vec t[m][m];
      vec_aliases_create2 (T, m, t); /* aliases to subsections */

      const real dr = rmax / nrad;

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          vec_rtab (HD, nrad, t_rad[i][j], dr, t[i][j]);

      vec_aliases_destroy2 (T, m, t);
    }
  }
#endif

  /*
    Most of these  Vecs will be a functional of  primary variable t or
    constant.  Earlier  versions supported  direct correlation c  as a
    primary variable.
  */
  local Vec c[m][m];
  vec_create2 (HD->da, m, c);

  /* Prepare intra-molecular correlations. The origin is at the corner
     as suitable for convolutions. Diagonal will be NULL: */
  local Vec w_fft[m][m];
  bgy3d_omega_fft_create (HD, m, solvent, w_fft); /* creates them */

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

  {
    local Vec t_fft[m][m];
    vec_create2 (HD->dc, m, t_fft); /* complex */

    local Vec c_fft[m][m];
    vec_create2 (HD->dc, m, c_fft); /* complex */

    /*
      Find a T such that dT as returned by iterate_t2 (&ctx, T, dT) is
      zero.   Cast is there  to silence  the mismatch  in the  type of
      first pointer argument: struct Ctx2* vs. void*:
    */
    {
      Ctx2 ctx =
        {
          .HD = HD,
          .m = m,
          .v_short = (void*) v_short,       /* in */
          .v_long_fft = (void*) v_long_fft, /* in */
          .y = (void*) c,                   /* work, c(t)*/
          .t_fft = (void*) t_fft,           /* work, fft(t) */
          .c_fft = (void*) c_fft,           /* work, fft(c(t)) */
          .w_fft = (void*) w_fft,           /* in */
        };

      bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_t2, T);
    }

    /* Free local stuff */
    vec_destroy2 (m, t_fft);
    vec_destroy2 (m, c_fft);
  }
  /*
    The approach with the primary variable  x == t won over time.  The
    direct correlation c(t) is  a dependent quantity! FIXME: we assume
    that the  value of c on  exit from the SNES  solver corresponds to
    the final t. If not, recompute it with compute_c() again.
  */

  Vec h[m][m];                  /* FIXME: not deallocated! */
  vec_create2 (HD->da, m, h);

  /* Compute h[][]: */
  {
    /* Now  that T  has  converged  we will  post-process  it using  the
       aliases t[][]: */
    local Vec t[m][m];
    vec_aliases_create2 (T, m, t); /* aliases to subsections */

    /* h = c + t: */
    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i; j++)
        VecWAXPY (h[i][j], 1.0, c[i][j], t[i][j]);

    vec_aliases_destroy2 (T, m, t);
  }
  vec_pack_destroy2 (&T);

  /* Chemical potential */
  {
    const real beta = HD->PD->beta;
    const real L3 = PD->L[0] * PD->L[1] * PD->L[2];

    /*
      FIXME: this is the only  place where one needs real-space rep of
      the  long-range  Coulomb.   So   far  it  is  computed  by  FFT.
      Alternative is to tablulate it  on the real-space grid using the
      analytic expression. Surprisingly, for a couple of tests we made
      (LJC,TIP3P)  the difference between  the two  approaches appears
      small.
    */
    local Vec cl[m][m];        /* real-space long-range correlation */
    vec_create2 (HD->da, m, cl);

    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i; j++)
        {
          /* Get real representation of long-range Coulomb potential */
          MatMultTranspose (HD->fft_mat, v_long_fft[i][j], cl[i][j]);

          /*
            Scale  Vec cl  to  get long-range  correlation (division  by
            volume is part of the inverse FFT):

              c  = -βv
               L      L
          */
          VecScale (cl[i][j], -beta/L3);
        }

    /* The function chempot() operates  on rectangular arrays, we pass
       an m x m square ones: */
    print_chempot (HD, m, m, h, c, cl);

    vec_destroy2 (m, cl);
  }

  /* No more used: */
  vec_destroy2 (m, c);
  vec_destroy2 (m, v_short);
  vec_destroy2 (m, v_long_fft);

  /*
    Compute  the  solvent susceptibility,  χ  =  ω  + ρh,  in  Fourier
    representation for  future use in  solute/solvent calculations. It
    appears to be more  convenient to operate with δχ = χ  - 1 so this
    is what is actually saved to disk.
  */
  {
    local Vec chi[m][m];
    vec_create2 (HD->da, m, chi);

    local Vec chi_fft[m][m];
    vec_create2 (HD->dc, m, chi_fft);

    const real dV = PD->h[0] * PD->h[1] * PD->h[2];

    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i; j++)
        {
          /*
              δχ = χ -  1 = (ω - 1) + ρh.

            The first ω-term is  only well representable on the k-grid
            and will be added later:
          */
          VecSet (chi[i][j], 0.0);
          VecAXPY (chi[i][j], PD->rho, h[i][j]);

          /* δχ(k) = FFT(δχ(r)) */
          MatMult (HD->fft_mat, chi[i][j], chi_fft[i][j]);
          VecScale (chi_fft[i][j], dV);

          /* Translate the  distribution to the grid  corner.  This is
             what one expects in convolution integrals: */
          bgy3d_vec_fft_trans (HD->dc, HD->PD->N, chi_fft[i][j]);

          /*
            Intra-molecular correlation  ω(k) is centered  at the grid
            corner, so add this only after translating chi_fft[][]:
          */
          if (w_fft[i][j])
            VecAXPY (chi_fft[i][j], 1.0, w_fft[i][j]);
          else
            /* ω - 1 == 0 identically, nothing to do! */
            assert (i == j);
        }

    /* The Fourier  rep is what is actually  read in solvent_kernel(),
       see below: */
    bgy3d_vec_save2 ("x%d%d-fft.bin", m, chi_fft);

    vec_destroy2 (m, chi);
    vec_destroy2 (m, chi_fft);
  }

  /* no more used: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)   /* sic! */
      vec_destroy (&w_fft[i][j]); /* Diagonal is NULL! */

  /* Delegated to the caller: vec_destroy (&h); */

  bgy3d_state_destroy (HD);

  /* g = 1 + h, store in Vec h for output: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecShift (h[i][j], 1.0);  /* FIXME: misnomer! */

  /* Dont forget to vec_destroy2() it! */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      g[i][j] = g[j][i] = h[i][j]; /* FIXME: misnomer! */

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

  /* Get solvent name or stay with the default: */
  bgy3d_getopt_string ("--solvent", name, sizeof name);

  /* Get the number of solvent sites and their parameters. Get it from
     the solute tables: */
  bgy3d_solute_get (name, &m, &solvent);

  /* Show solvent parameters: */
  bgy3d_sites_show ("Solvent", m, solvent);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", name);

  Vec g[m][m];
  hnc3d_solvent_solve (PD, m, solvent, g);

  vec_destroy2 (m, g);

  return PETSC_NULL;
}


/*
  In k-space compute either

    t   = (χ   - 1) * c
     vu     vv         vu

  or

    t   = ρ c   *  h
     vu      vv     vu

  Compare to Eqs.   (7) and (8) in the  header comments, respectively.
  Given that the solute index  "u" is redundant, both expressions have
  the structure of the matrix-vector product of a (fixed) solvent pair
  quantitiy  and a  (variable) "vector"  of site  distributions around
  solute    impurity    to    give    another   "vector"    of    site
  distributions. FIXME: the names of the local variables correspond to
  Eq. (8)  for historical reasons ---  the variant of  Beglov and Roux
  was  there first.   Also note  that the  order of  multiplication χc
  assumes that χ is symmetric which is only the case if number density
  of all solvent sites is the same.

  Here ρ  is a scalar overall  factor equal to solvent  density in one
  case or just 1 in the other case (only if the convolution theorem is
  factorless).  The caller may (and does) abuse this factor to further
  scale the result.

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
  Vec (*chi_fft)[m] = (void*) ctx->c_fft; /* [m][m], complex, in */
  Vec *c = (void*) ctx->y;                /* [m], real, work */
  Vec *c_fft = (void*) ctx->h_fft;        /* [m], complex, work */

  const ProblemData *PD = ctx->HD->PD;
  // const real rho = PD->rho;
  const real beta = PD->beta;
  const real N3 = PD->N[0] * PD->N[1] * PD->N[2];

  /* Establish aliases to the subsections of the long Vec T and dT: */
  local Vec t[m], dt[m];
  vec_aliases_create1 (T, m, t);
  vec_aliases_create1 (dT, m, dt);

  /*
    The  new  candidate for  the  direct  uv-correlation  c using  the
    closure relation:
  */
  for (int i = 0; i < m; i++)
    {
      compute_c (PD->closure, beta, ctx->v[i], t[i], c[i]);

      /* fft(c).  Here  c is the  3d unknown direct  uv-correlation of
         the solvent sites and the solute species as the whole: */
      MatMult (ctx->HD->fft_mat, c[i], c_fft[i]);
    }

  /*
    Now

      t   = (χ - 1) * c
       vu     vv       vu

    Here χ  is the constant solvent  susceptibility.  The "convolution
    star" corresponds  to a matrix multiplication in  the k-space: for
    each solvent site sum over  solvent sites.  Note that the input in
    chi_fft[][] corresponds to χ - 1.

    Let the overall  scale include a factor for  forward & inverse FFT
    right away  --- check how  we (i) conveniently forgot  to multiply
    c_fft[] with grid weight h³  after forward FFT and (ii) the result
    of inverse FFT is also not divided by L³ as everywhere else:
  */
  compute_t1 (m, 1.0 / N3, chi_fft, c_fft, ctx->t_fft);

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


/* Reads  χ -  1 into  chi_fft[][] as  previousely written  by solvent
   solver.  See hnc3d_solvent_solve() above. */
static void solvent_kernel (State *HD,
                            int m, const Site solvent[m], /* in */
                            Vec chi_fft[m][m])            /* out */
{
  if (bgy3d_getopt_test ("--from-radial-g2")) /* FIXME: better name? */
    {
      /*
        Load  radial   data  from  text   file.   FIXME:  representing
        long-range on a real-space grid is intrinsically broken!
      */
      const real dV = HD->PD->h[0] * HD->PD->h[1] * HD->PD->h[2];

      local Vec chi[m][m];
      vec_create2 (HD->da, m, chi);

      bgy3d_vec_read_radial2 (HD->da, HD->PD, "x%d%d.txt", m, chi);

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          {
            MatMult (HD->fft_mat, chi[i][j], chi_fft[i][j]);
            VecScale (chi_fft[i][j], dV);

            /* Translate the distribution to  the grid corner. This is
               what one expects in convolution integrals. */
            bgy3d_vec_fft_trans (HD->dc, HD->PD->N, chi_fft[i][j]);
          }
      vec_destroy2 (m, chi);
    }
  else
    bgy3d_vec_read2 ("x%d%d-fft.bin", m, chi_fft); /* ready for use as is */

#ifdef WITH_FORTRAN
  if (false)
    {
      /*
        Problem data for  use in 1d-RISM code. If  you do not upscale,
        the rmax =  max (PD->L) / 2 will be  too low for interpolation
        on the r-grid. Though here only the k-grid is of interest.
      */
      ProblemData pd = rism_upscale (HD->PD);

      const int nrad = rism_nrad (&pd);
      const real rmax = rism_rmax (&pd);

      real x_fft[m][m][nrad];

      /* 1d-RISM  calculation. Dont  need  indirect correlation,  need
         χ(k): */
      rism_solvent (&pd, m, solvent, NULL, x_fft);

      /* The solute/solvent code expects χ - 1: */
      for (int i = 0; i < m; i++)
        for (int k = 0; k < nrad; k++)
          x_fft[i][i][k] -= 1;

      /* dr * dk = π / nrad where dr = rmax / nrad */
      const real dk = M_PI / rmax;

      /* FIXME: assuming symmetric χ - 1 with aliasing: */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          vec_ktab (HD, nrad, x_fft[i][j], dk, chi_fft[i][j]);
    }
#endif
}


/*
  Solving  HNC3d equations.   Solvent  susceptibility χ  -  1 of  pure
  solvent  appears  as a  fixed  input  here.  A primary  variable  is
  indirect correlation t related to h by one of the closure relations.

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
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         Vec g[m],
                         Context **medium)
{
  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  PetscPrintf (PETSC_COMM_WORLD, "(iterations for γ)\n");

  State *HD = bgy3d_state_make (PD); /* FIXME: rm unused fields */

  /*
    This will be a functional h(t) of primary variable:
  */
  Vec h[m];                     /* FIXME: not deallocated */
  vec_create1 (HD->da, m, h);

  local Vec v[m];
  vec_create1 (HD->da, m, v); /* solute-solvent interaction */

  /*
    Get  solute-solvent interaction.  Fill  v[i] with  the short-range
    potential acting on solvent site "i". FIXME: asymptotic Coulomb is
    assumed to be short-range and is added to to total potential:
  */

  /* Asymptotic Coulomb field of the (often neutral) solute: */
  local Vec uc = vec_create (HD->da);
  local Vec uc_rho = vec_create (HD->da); /* FIXME: discarded */

  bgy3d_solute_field (HD, m, solvent, n, solute,
                      v, uc,  /* out */
                      uc_rho, /* smeared core density, discarded */
                      density);  /* void (*density)(...) */

  /*
    If   the  solute  is   neutral,  asymptotic   electrostatics  is
    short-range  (if  you  dare   to  call  dipole  field  that,  of
    course). Add it to the rest:
  */
  for (int i = 0; i < m; i++)
    VecAXPY (v[i], solvent[i].charge, uc);

  /*
    For primary variable t there are  two ways to access the data: via
    the long Vec  T and m shorter Vec t[m]  aliased to the subsections
    of the longer one.
  */
  local Vec T = vec_pack_create1 (HD->da, m); /* long Vec */

  /* Zero as intial guess for t  is (almost?) the same as exp(-βv) - 1
     initial guess for h: */
  VecSet (T, 0.0);

  {
    local Vec h_fft[m];
    vec_create1 (HD->dc, m, h_fft); /* complex */

    local Vec t_fft[m];
    vec_create1 (HD->dc, m, t_fft); /* complex */

    /* This should be the only pair quantity, χ - 1: */
    local Vec chi_fft[m][m];
    vec_create2 (HD->dc, m, chi_fft);        /* complex */

    /*
      Get the solvent-solvent  susceptibility (offset by one).  FIXME:
      we get the solvent description as an argument, but read the file
      to  get the rest.   There is  no guarantee  the two  sources are
      consistent.
    */
    solvent_kernel (HD, m, solvent, chi_fft);

    /*
      Find T  such that dt  as returned by  iterate_t1 (HD, T,  dT) is
      zero. Cast is there to silence the mismatch in the type of first
      pointer argument: Ctx1* vs. void*:
    */
    {
      /* Work area for iterate_t1(): */
      Ctx1 ctx =
        {
          .HD = HD,
          .m = m,
          .v = v,
          .y = h,                   /* h(t), not t(h) */
          .c_fft = (void*) chi_fft, /* pair quantitity */
          .h_fft = h_fft,
          .t_fft = t_fft,
        };

      bgy3d_snes_default (PD, &ctx, (VectorFunc) iterate_t1, T);
    }
    vec_destroy1 (m, h_fft);
    vec_destroy1 (m, t_fft);

    /* This should have been the only pair quantity: */
    vec_destroy2 (m, chi_fft);
  }

  /*
    Now that the solution T has  converged, use it to compute the rest
    of the  quantities.  (Re)evaluate  c(t) and h(t)  without assuming
    that e.g. Vec h already  has any meaningful value.  We need direct
    correlation c to compute chemical potential:
  */
  local Vec c[m];
  vec_create1 (HD->da, m, c);

  {
    local Vec t[m];
    vec_aliases_create1 (T, m, t); /* aliases to subsections */

    /* Closure relation, c = c(t): */
    for (int i = 0; i < m; i++)
      compute_c (PD->closure, PD->beta, v[i], t[i], c[i]);

    /* h = c + t */
    for (int i = 0; i < m; i++)
      VecWAXPY (h[i], 1.0, c[i], t[i]);

    vec_aliases_destroy1 (T, m, t);
  }
  vec_pack_destroy1 (&T);       /* not vec_destroy()! */
  vec_destroy1 (m, v);

  /* Excess chemical potential: */
  {
    /*
      So  far  the  solute/solvent  code does  not  handle  long-range
      potential of the solute (assuming  it is neutral).  The code for
      the chemical potential expects that. Supply zeroes:
    */
    local Vec cl[m];
    vec_create1 (HD->da, m, cl);

    for (int i = 0; i < m; i++)
      VecSet (cl[i], 0.0);

    /*
      In  3d  models  the  solute  is  effectively  a  single  (albeit
      non-spherical) site. Treat the arrays  h[m], c[m] and cl[m] as 1
      x m arrays here:
    */
    print_chempot (HD, 1, m, (void*) h, (void*) c, (void*) cl);

    vec_destroy1 (m, cl);
  }
  /* No more used: */
  vec_destroy1 (m, c);

  /* g := h + 1 */
  for (int i = 0; i < m; i++)
    VecShift (h[i], 1.0);       /* FIXME: misnomer! */

  /* return the Context **ret to caller for integration over electron
   * potential of solvent */
  Context *ret = info (HD, m, solvent, n, solute, h, uc, uc_rho);
  if (medium)
    *medium = ret;
  else
    bgy3d_pot_destroy (ret);

  /* keep uc and uc_rho until being used in info() */
  vec_destroy (&uc);
  vec_destroy (&uc_rho);

  /* free stuff */
  /* Delegated to the caller: vec_destroy (&h) */

  bgy3d_state_destroy (HD);

  /* Caller is supposed to destroy it! */
  for (int i = 0; i < m; i++)
    g[i] = h[i];                /* FIXME: misnomer! */

  bgy3d_vec_save1 ("g%d.bin", m, g);
}


/*
  Solving for  h of HNC equation  with a default  non-linear solver as
  specified by the --snes-solver option.  The solvent susceptibility χ
  - 1 is fixed and appears as an input here:
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

  /* Get solvent name or stay with the default: */
  bgy3d_getopt_string ("--solvent", name, sizeof name);

  /* Get the number of solvent sites and their parameters. Get it from
     the solute tables: */
  bgy3d_solute_get (name, &m, &solvent);

  /* Show solvent parameters: */
  bgy3d_sites_show ("Solvent", m, solvent);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solvent is %s.\n", name);

  /* Get solute name or stay with the default: */
  bgy3d_getopt_string ("--solute", name, sizeof name);

  /* Get the solute from the tables: */
  bgy3d_solute_get (name, &n, &solute);

  /* Show solute parameters: */
  bgy3d_sites_show ("Solute", n, solute);

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solute is %s.\n", name);

  Vec g[m];
  hnc3d_solute_solve (PD, m, solvent, n, solute,
                      NULL,  /* no electron density */
                      g,
                      NULL); /* no medium returned */

  vec_destroy1 (m, g);
  return PETSC_NULL;
}
