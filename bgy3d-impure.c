/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2OS.c,v 1.20 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

/*
  According to Page. 115 - 116:

    g(x) = g (x) exp[-u(x)]                                        (1)
            0

  with

    g (x) = exp[-β (V  (x) + V  (x) + V  (x))]                     (2)
     0               LJ       CS       CL

  Here  and below  LJ stays  for Lennard-Jones,  CS stays  for Coulomb
  short, and  CL stays for Coulomb  long.  Then g(x)  can be rewritten
  as:

           ~
    g(x) = g  (x) exp[-β V  (x) - u(x)]
             0            CL

  with

    ~
    g  (x) = exp[-β (V  (x) + V  (x))]
      0               LJ       CS

  then BGY equation written as:

     ~
    Δu = K(g) + β ΔV
                    CL

  with
           ~                             ~          ~
    g(x) = g (x) exp[-β V  (x) - u(x)] = g (x) exp[-u(x)]
            0            CL               0

                  ~
  the solution of u can be represented by a difference of two functions:

    ~   -    *
    u = u - u

  while:

     -                       -
    Δu = K(g) + β ΔV   in Ω, u(∂Ω) = f,
                    CL

      *           *
    Δu = 0 in Ω, u (∂Ω) = f,

  so after solving:

    Δu = K(g)

  we get:

    -
    u = u + β V
               CL

  Just  like  in  the  last  equation  the  original  code  added  the
  long-range  Coulomb term  with the  (back then  erroneously missing)
  inverse temperature  factor to  the final PMF,  u(x), at the  end of
  iteration  after computing the  convolution term  and all  the magic
  intra-molecular  contributions. By now  it is  treated more  like in
  Eqs.   (1)  and  (2) where  it  is  a  part  of zero  density  limit
  potential. Instead the the initial  guess for PMF u(x) is defined to
  (almost) negative of the Coulomb long: - [(ε - 1) / ε] uc(x).  There
  is a kind of arbitrarness  whether to assign constant terms to g₀(x)
  or to u(x).

  The Coulomb  potential of  the solute common  for all  solvent sites
  that  differ only  by their  charges is  calculated  beforehand. The
  site-specific potentials  can be obtained by  multiplying the common
  long-range Coulomb potential (stored in  Vec uc) by the solvent site
  charge.

  There is  another story about short- and  long-range Coulomb. Namely
  that of the two such terms in the *solvent* kernel K[]. The two were
  intrinsically related in the pure  code where solute and the solvent
  are the same species. FIXME: elaborate this point.
*/


#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* needs Site */
#include "bgy3d-force.h"        /* bgy3d_force() */
#include "bgy3d-getopt.h"
#include "bgy3d-vec.h"
#include "bgy3d-pure.h"         /* bgy3d_nssa_intra_log() */
#include "bgy3d-poisson.h"      /* bgy3d_poisson() */
#include "bgy3d-dirichlet.h"    /* Laplace staff */
#include "bgy3d-potential.h"    /* Context, etc. */
#include "bgy3d-fftw.h"         /* bgy3d_fft_interp() */
#include "bgy3d-snes.h"         /* bgy3d_snes_newton() */
#include "bgy3d-impure.h"

#include <stdbool.h>            /* bool */
#include <complex.h>            /* after fftw.h */

/* FIXME: bgy3d-solvents.h pollutes the namespace: */
#undef sH
#undef eH
#undef qH
#undef sO
#undef eO
#undef qO
#undef G
#undef r_HO

/* Side effects: uses BHD->fg2_fft[3] as work arrays. */
static void  pair (State *BHD,
                   const Site a, const Site b, /* struct by value? */
                   const Vec g2,
                   Vec f_short[3], Vec f_long[3], /* work arrays */
                   Vec fs_g2_fft[3],              /* intent(out) */
                   Vec fl_g2_fft[3],              /* intent(out) */
                   Vec u2, Vec u2_fft)            /* intent(out) */
{
  /* Compute  forces  and long-range  Coulomb  (in  both  reps) for  a
     pair: */
  bgy3d_force (BHD, a, b, f_short, f_long, NULL, NULL, u2, u2_fft,
               1.0, 1.0);

  /*
    Compute weighted forces FFT(F *  g^2). Here star is a product, not
    a convolution.

           (2)                 (2)          (2)
      F * g   = (F   + F  ) * g    + (F  * g   - F  ) + F
                  LJ    CS             CL         CL     CL

    (LJ is Lennard-Jones, CS is Coulomb short, CL is Coulomb long).

    See (5.101) and (5.102)  in Jager's thesis.  The k-representations
    of  the  first  two  terms  will be  stored  in  fs_g2_fft[3]  and
    fl_g2_fft[3], respectively.  Note that the latter, fl_g2_fft[], is
    a long-range  force weighted by  g2 - 1  and not g2.   The missing
    FFT(F_CL)  has  been  calculated  by  ComputeFFTfromCoulomb()  but
    discarded.  The code needs at least one work vector, use this:
  */
  local Vec work = vec_pop (BHD->da);

  FOR_DIM
    {
      /* First (F_LJ + F_CS) * g2: */
      VecPointwiseMult (work, g2, f_short[dim]);

      /* Next FFT((F_LJ + F_CS) * g2): */
      MatMult (BHD->fft_mat, work, fs_g2_fft[dim]);

      /* Now Coulomb long, F_CL * g2: */
      VecPointwiseMult (work, g2, f_long[dim]);

      /* Next F_CL * g2 - F_CL: */
      VecAXPY (work, -1.0, f_long[dim]);

      /* Finally FFT(F_CL * g2 - F_CL): */
      MatMult (BHD->fft_mat, work, fl_g2_fft[dim]);
    }

  vec_push (BHD->da, &work);
}

/*
  This is supposed to compute  the k-grid representation of the matrix
  (linear operator) that relates unknown u to the RHS g:

    Δu = - β ρ K * g, K = div (F g2)

  with "*"  being the convolution  that in Fourier rep  degenerates to
  pointwise  product,  g2  and  g being  the  (fixed)  solvent-solvent
  site-site  distribution  function  and  (the unknown)  solvent  site
  distribution around the solute.

  The kernel dfg(kx,  ky, kz) includes the Poisson  factor as in -K/k2
  but does not include the factor βρ:

    u(kx, ky, kz) = β ρ dfg(kx, ky, kz) g(kx, ky, kz)

  The three  fg[3] arrays  are intent(in).  The  dfg array  for kernel
  -K/k2 is intent(out).

  The  weighted force  fg[3] is  missing the  long range  Coulomb part
  F_CL,  see comments  in pair()  above. Note  that divergence  of the
  gradient is the  laplacian and F_CL is the  negative of the gradient
  of the long-range  Coulomb interaction: div (F_CL) =  - Δu2. That is
  why the  function accepts additional  Coulomb potential in  Vec coul
  and adds that (almost) unmodified to the convolution kernel.

  No known side effects.
*/
static void kernel (const DA dc,
                    const ProblemData *PD,
                    Vec fg[3],  /* complex, intent(in) */
                    Vec coul,   /* complex, intent(in) */
                    Vec dfg)    /* complex, intent(out) */
{
  int x[3], n[3], i[3];

  const int *N = PD->N;         /* [3] */
  const real *L = PD->L;        /* [3] */

  /* FIXME: rectangular box? */
  assert (L[0] == L[1]);
  assert (L[0] == L[2]);

  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L3 = L[0] * L[1] * L[2];
  const real fac = L[0] / (2.0 * M_PI); /* BHD->f ist nur grad U, nicht F=-grad U  */
  const real scale = fac / L3;

  /* Get local portion of the grid */
  DMDAGetCorners (dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* Loop over local portion of grid: */
  complex ***fg_[3], ***dfg_, ***coul_;
  DMDAVecGetArray (dc, dfg, &dfg_);
  if (coul)
    DMDAVecGetArray (dc, coul, &coul_);
  FOR_DIM
    DMDAVecGetArray (dc, fg[dim], &fg_[dim]);

  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          real k2 = SQR(ic[2]) + SQR(ic[1]) + SQR(ic[0]);

          real k_fac;
          if (k2 > 0)
              k_fac = h3 * h3 * scale / k2;
          else
              k_fac = 0.0;

          /*
            Compute the  (Fourier transform  of) of the  divergence of
            the "weighted  force" vector div  (F g) which serves  as a
            convolution  kernel  in BGY3dM  equations.   Note the  the
            derivative in momentum space involves a factor -I:

                div = -I dot(k, (F g))

            The  additional  factor  -k^(-2) effectively  included  in
            k_fac recovers the  Fourier transform of the corresponding
            Poisson solution.
           */

          /*
            Complex arithmetics here.  "I"  is a macro expanding to an
            imaginary unit:
          */
          complex sum = 0.0;
          for (int p = 0; p < 3; p++)
            sum += ic[p] * fg_[p][i[2]][i[1]][i[0]];

          dfg_[i[2]][i[1]][i[0]] = k_fac * (-I * sum);

          /*
            Origin of this occasional addition is the custom treatment
            of  the  long-range Coulomb  interaction  that avoids  the
            round  trip:  potential  ->  gradients  ->  divergence  ->
            inverse laplacian.  Note that F_CL force is missing in the
            long-range weighted  forces which include only  F_CL (g2 -
            1).  See also description of the function.

            The   difference    between   Compute_H2O_interS_C()   and
            Compute_H2O_interS() of the  original code reduces to this
            addendum.  Supply a NULL pointer for Vec coul to omit it.

            Long range Coulomb part (right one):
          */
          if (coul)
            dfg_[i[2]][i[1]][i[0]] += (h3 / L3) * coul_[i[2]][i[1]][i[0]];
        }
  DMDAVecRestoreArray (dc, dfg, &dfg_);
  if (coul)
    DMDAVecRestoreArray (dc, coul, &coul_);
  FOR_DIM
    DMDAVecRestoreArray (dc, fg[dim], &fg_[dim]);

  /*
    Translate  the kernel  so  that the  real-space representation  is
    centered at the  grid corner (0, 0, 0) and not  at the grid center
    like other  grid representations. This  is what is assumed  in the
    convolution integrals:
  */
  bgy3d_vec_fft_trans (dc, N, dfg);
}


/*
  This applies the  kernel compured by kernel() to FFT  of g to obtain
  an increment to du. Such a product in the k-space

    du  := du  + scale * K  * g
      k      k            k    k

  corresponds to convolution in the real space.

  The kernel was derived from  the Fourier transform of the divergence
  of  the "weighted  force" vector  div (F  g). The  additional factor
  -k^(-2)  effectively included  in  the kernel  recovers the  Fourier
  transform of the corresponding  Poisson solution. The factor scale =
  βρ is not included, on the other hand. See kernel() for details.

  The convolution  kernel in BGY3dM  serves role similar to  that pure
  solvent direct correlation c in HNC3d equations.

  Complex  Vec du  is intent(inout).   Dont forget  to clear  it early
  enough.
*/
static void apply (Vec ker_fft, Vec g_fft, const real scale, /* in */
                   Vec du_fft)  /* inout */
{
  /* FMA stays for "fused multiply-add" */
  complex pure fma (complex du, complex c, complex g)
  {
    return du + scale * (c * g);
  }

  vec_fft_map3 (du_fft, fma, du_fft, ker_fft, g_fft); /* aliasing! */
}

/*
  Given  the  solvent  description  and  pair  distribution  functions
  compute the multiplicative  complex momentum-space representation of
  kernel to be applied every iteration.

  Side  effects:  by  way  of  pair() uses  BHD->fg2_fft[3]  as  work
  arrays.
*/
static void solvent_kernel (State *BHD, int m, const Site solvent[m],
                            Vec g2[m][m],         /* in */
                            Vec kernel_fft[m][m]) /* out */
{
  /* Real  work  vectors,  re-used   for  all  m(m+1)/2  solvent  site
     pairs. */
  local Vec u2 = vec_create (BHD->da);
  local Vec fs[3], fl[3];
  vec_create1 (BHD->da, 3, fs); /* 3-vector */
  vec_create1 (BHD->da, 3, fl); /* 3-vector */

  /* Complex work vectors, re-used for all pairs: */
  local Vec kl_fft = vec_create (BHD->dc);
  local Vec u2_fft = vec_create (BHD->dc);
  local Vec fs_g2_fft[3], fl_g2_fft[3];
  vec_create1 (BHD->dc, 3, fs_g2_fft); /* 3-vector */
  vec_create1 (BHD->dc, 3, fl_g2_fft); /* 3-vector */

  /* Over all pairs: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* Does real work: */
        pair (BHD, solvent[i], solvent[j], g2[i][j],
              fs, fl, /* work vectors*/
              fs_g2_fft, fl_g2_fft, u2, u2_fft);

        /* Kernel term originating from the the short-range (weighted)
           forces: */
        kernel (BHD->dc, BHD->PD, fs_g2_fft, NULL, kernel_fft[i][j]);

        /* This used to be  the long-range contribution, we compute it
           separately, but later add to the total kernel: */
        kernel (BHD->dc, BHD->PD, fl_g2_fft, u2_fft, kl_fft);

        /*
          FIXME:  what is  the point  to split  the kernel  in two
          pieces?  Redefine S := S + L and forget about L:

          ker = ker_l + ker_s
        */
        VecAXPY (kernel_fft[i][j], 1.0, kl_fft);
      }

  /* Clean up and exit: */
  vec_destroy (&u2);
  vec_destroy (&u2_fft);
  vec_destroy (&kl_fft);
  vec_destroy1 (3, fs);
  vec_destroy1 (3, fl);
  vec_destroy1 (3, fs_g2_fft);
  vec_destroy1 (3, fl_g2_fft);
}


/* g := exp (-u) */
static void mexp (Vec g, Vec u)
{
  real pure f (real u)
  {
    return exp (-u);
  }
  vec_map1 (g, f, u);
}

static void iterate (State *BHD,
                     int m,
                     const Site solvent[m],
                     const real rhos[m],
                     Vec kernel_fft[m][m], /* in */
                     Vec omega_fft[m][m],  /* in */
                     Vec u0[m],            /* in */
                     Vec uc,               /* in */
                     Vec u[m],             /* in */
                     Vec x_lapl[m],        /* inout */
                     Vec g[m],             /* out */
                     Vec g_fft[m],         /* work */
                     Vec du_acc_fft,       /* work */
                     Vec work,             /* work */
                     Vec du[m])            /* out */
{
  /* Some functions,  such as bgy3d_nssa_intra_log()  use preallocated
     complex Vecs in State BHD for work. */

  const ProblemData *PD = BHD->PD;
  const real beta = PD->beta;

  /* Apply  boundary  conditions  (or  other constraints)  and  update
     g[]. Intent(out) Vec du[] can be used as work arrays: */
  for (int i = 0; i < m; i++)
    {
      /*
        Set  accumulator to  the  sum of  short-range  field u0[]  and
        long-range Coulomb field  uc scaled by the site  charge to the
        accumulator.   Beware that the  original code  is erroneousely
        missing the inverse temperature beta in this expression:
      */
      VecWAXPY (du[i], beta * solvent[i].charge, uc, u0[i]);

      VecAXPY (du[i], 1.0, u[i]);

      /*
        See  pp.  116-177  in  thesis: boundary  conditions (5.107)  -
        (5.110): first impose  boundary condition then solve laplacian
        equation  and  subtract from  the  field.   State  BHD is  not
        modified by these calls.

        Vec du[] and x_lapl[]  are intent(inout) here. Try to preserve
        the values of  x_lapl[] across iterations to save  time in the
        iterative solver.
      */
      bgy3d_impose_laplace_boundary (BHD, du[i], x_lapl[i]);
    }

  /* g := exp (-u): */
  for (int i = 0; i < m; i++)
    mexp (g[i], du[i]);

  /* Compute FFT of g[] for all sites: */
  for (int i = 0; i < m; i++)
    MatMult (BHD->fft_mat, g[i], g_fft[i]);

  /* for H, O in that order ... */
  for (int i = 0; i < m; i++)
    {
      /*
        ... sum  over H, O in  that order.  LJ, short-  and long range
        Coulomb, and a so called  strange addition is accounted for in
        the kernel.  First clear accumulator:
      */
      VecSet (du_acc_fft, 0.0);

      for (int j = 0; j < m; j++) /* This increments the accumulator: */
        apply (kernel_fft[i][j], g_fft[j], beta * rhos[j], du_acc_fft);

      /*
        Compute  IFFT  of  du_acc_fft  for the  current  site.   Other
        contributions are added to the real space du[] below:
      */
      MatMultTranspose (BHD->fft_mat, du_acc_fft, du[i]);

      /*
        In  the following the  sum is  over all  sites j  /= i.  It is
        supposed to account  for intra-molecular site-site correlation
        of   the  solvent  sites   due  to   the  rigid   bonds.   See
        e.g. Eq. (4.114), Jager Diss:
      */
      for (int j = 0; j < m; j++)
        {
          if (j == i) continue;
          /*
            Compute  intramolecular contribution  due to  j /=  i.  So
            called  log-term, see  e.g.  Eqs.   (4.107),  (4.109), and
            (4.114).  Implies that the  two sites  belong to  the same
            species.  This  body is  executed exactly once  for 2-site
            models.  Vec   work  is  intent  (out)  to   be  added  to
            accumulator later:
          */
          bgy3d_nssa_intra_log (BHD, g_fft[i], omega_fft[i][j], g[j],
                                work); /* out */

          /* Add the contribution of site j /= i to du[i]: */
          VecAXPY (du[i], 1.0, work);

          /*
            FIXME:  the  code  for   pure  solvent  adds  a  so-called
            "numerically  challenging"  term  due  to  intra-molecular
            correlations too.  There seem to be  no corresponding term
            in this code.
          */
        }

      /* By convention,  we return the  increment, not the new  u[] so
         that at convergence the result vanishes: */
      VecAXPY (du[i], -1.0, u[i]);
    } /* over sites i */
}


typedef struct Ctx
{
  State *BHD;
  int m;
  const Site *solvent;          /* [m] */
  real *rhos;                   /* [m] */
  Vec *kernel_fft;              /* [m][m] */
  Vec *omega_fft;               /* [m][m] */
  Vec *u0;                      /* [m] */
  Vec uc;                       /*  */
  Vec *x_lapl;                  /* [m] */
  Vec *g;                       /* [m] */
  Vec *g_fft;                   /* [m] */
  Vec du_acc_fft;               /* work */
  Vec work;                     /* work */
} Ctx;


static void iterate_u (Ctx *s, Vec U, Vec dU)
{
  const int m = s->m;           /* number of solvent sites */

  /* Compiler expects us to pass [m][m] arrays to iterate(): */
  Vec (*kernel_fft)[m] = (void*) s->kernel_fft;
  Vec (*omega_fft)[m] = (void*) s->omega_fft;

  /* Establish aliases to the subsections of the long Vec U and dU: */
  local Vec u[m], du[m];
  vec_aliases_create1 (U, m, u);
  vec_aliases_create1 (dU, m, du);

  /* Iterate u[] -> du[]: */
  iterate (s->BHD,              /*  1 in */
           s->m,                /*  2 in */
           s->solvent,          /*  3 in */
           s->rhos,             /*  4 in */
           kernel_fft,          /*  5 in */
           omega_fft,           /*  6 in */
           s->u0,               /*  7 in */
           s->uc,               /*  8 in */
           u,                   /*  9 in <- here */
           s->x_lapl,           /* 10 inout */
           s->g,                /* 11 out */
           s->g_fft,            /* 12 work */
           s->du_acc_fft,       /* 13 work */
           s->work,             /* 14 work */
           du);                 /* 15 out <- here */

  /*
    This  destroys the  aliases,  but  does not  free  the memory,  of
    course. The actuall data is owned by Vec U and Vec dU. From now on
    one may access U and dU directly again.
  */
  vec_aliases_destroy1 (U, m, u);
  vec_aliases_destroy1 (dU, m, du);
}


/*
 * FIXME: when  passing in a  NULL pointer, this destroy  function can
 * still be executed without error.
 *
 * From time to time we have to pass an artificial value of 'restart',
 * e.g.  hnc3d_solute_solve() neither  receives 'restart'  pointer nor
 * returns it,  but when  calling from PG,  an argument  for 'restart'
 * from  the scheme interface  'bgy3d_solute_hook' is  still expected.
 * In this case, we could (only)  set restart as '0' and since it will
 * not be refreshed during  calling hnc3d_solute_solve(), we will have
 * a   NULL  pointer  when   this  destroy   function  is   called  by
 * bgy3d_finalize() in PG.
 *
 * Due to the type cast  in bgy3d-guile.c, 'restart' could only be set
 * as integer and nothing else.
 */
void bgy3d_restart_destroy (Restart *restart)
{
  /* Restart info  is just a (long)  Vec that happens to  fit into (or
     *rather is*) a valid pointer or NULL: */
  local Vec U = (Vec) restart;
  if (U)
    vec_pack_destroy1 (&U);
}


/*
  This function  solves the the  BGY3dM equation for a  m-site solvent
  and  an arbitrary  solute.  Solvent  properties in  the form  of the
  array of  sites, convolution kernel  and intra-molecular correlation
  functions  should  be  supplied  by the  caller.  The  pre-allocated
  vectors

  Vec g[m], intent(out)

  provided   by  the  caller   are  filled   with  the   solvent  site
  distributions.

  Context **v, inent(out);

  Is set to  an iterator over the solvent  potential field. The caller
  is responsible  for calling bgy3d_pot_destroy() when it  is not more
  needed.
*/
static void solute_solve (State *BHD,
                          int m, const Site solvent[m], /* in */
                          Vec kernel_fft[m][m],         /* in */
                          Vec omega_fft[m][m],          /* in */
                          int n, const Site solute[n],  /* in */
                          void (*density)(int k, const real x[k][3], real rho[k]),
                          Vec g[m],          /* out */
                          Context **medium,  /* out, optional */
                          Restart **restart) /* inout, optional */
{
  const ProblemData *PD = BHD->PD;

  /* Site specific  density.  Computed as a solvent  density rho times
     number of sites of that type in a solvent: */
  real rhos[m];
  for (int i = 0; i < m; i++)
    rhos[i] = PD->rho;          /* 2 * PD->rho; */

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Inverse temperature: */
  const real beta = PD->beta;

#ifdef L_BOUNDARY
  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem across iterations. See call to KSPSetInitialGuessNonzero()
    in bgy3d-poisson.c:
  */
  local Vec x_lapl[m];          /* real */
  vec_create1 (BHD->da, m, x_lapl);
  for (int i = 0; i < m; i++)
    VecSet (x_lapl[i], 0.0);
#endif

  /*
    These complex  vectors will  hold FFT of  the current  g. Allocate
    enough to  hold a  local portion  of the grid  and free  after the
    loop.
  */
  local Vec g_fft[m];
  vec_create1 (BHD->dc, m, g_fft); /* complex */

  /*
    There is no  point to transform each contribution  computed in the
    momentum space back, accumulate them on the k-grid in this complex
    Vec:
  */
  local Vec du_acc_fft = vec_create (BHD->dc); /* complex */

  local Vec work = vec_create (BHD->da);

  /* Coulomb long, common for all sites: */
  local Vec uc = vec_create (BHD->da);

  /* Electron density for integration: */
  local Vec uc_rho = vec_create (BHD->da);

  /*
    Later u0  = beta *  (VM_LJ + VM_coulomb_short), which  is -log(g0)
    actually.  See: (5.106)  and (5.108) in Jager's thesis.  It is not
    filled with data yet, I assume.
  */
  local Vec u0[m];                    /* real */
  vec_create1 (BHD->da, m, u0); /* real */

  /*
    For primary  variable there are two  ways to access  the data: via
    the long Vec and m shorter  Vecs aliased to the subsections of the
    longer one.
  */
  local Vec U = vec_pack_create1 (BHD->da, m); /* long Vec */

  for (real damp = damp_start; damp <= 1.0; damp += 0.1)
    {
      /*
        Fill  u0[0], u0[1]  (see  the definition  above)  and uc  with
        VM_Coulomb_long.  No  other fields of the  struct State except
        those passed explicitly are modified:
      */
      bgy3d_solute_field (BHD, m, solvent, n, solute,
                          u0, uc, uc_rho, /* out */
                          density);       /* void (*density)(...) */

      /* Scale solute-solvent interactions: */
      PetscPrintf (PETSC_COMM_WORLD,
                   "Scaling solute-solvent interactions by %f\n", damp);
      /*
        Historically  short-range potential is  scaled by  the inverse
        temperature and the code operates  with u(x) = βv(x) insead of
        potential  itself   at  many  places.    FIXME:  actually  all
        potentials should appear only in this combination!
      */
      VecScale (uc, damp);
      VecScale (uc_rho, damp);
      for (int i = 0; i < m; i++)
        VecScale (u0[i], beta * damp);

      /*
        Set initial  guess, either  as supplied by  the caller,  by an
        estimate or  by reading  from file. At  the end of  the "dump"
        loop du[] is  (eventually) written to disk, so  that next time
        we can read an updated version.

        The  unscreened  Coulomb  field  of  the  solute  may  be  too
        attractive in  some regions close  to the excluded  volume, so
        that initial  distribution g0(x) has  unphysically large peaks
        in  that region.  In  such cases  (e.g.  CS2/QM)  the solution
        immediately  diverges.   The  original  code handled  this  by
        omitting the Coulomb long from g0(x) and added that to u(x) in
        the body of iterations.  We  here start with the initial guess
        u(x) = - [(ε - 1) / ε] * q * uc(x) with some large ε instead.
      */
      if (restart && *restart)
        {
          /*
            If the argument is present  and valid, this is the data to
            be used for resuming iterations.  So far this data is just
            a ref to a long Vec that happens to fit into a pointer:
          */
          Vec U_old = (Vec) *restart;

          /* Initialize long  Vec U by  copying restart data  from the
             last run: */
          VecCopy (U_old, U);
          bgy3d_restart_destroy (*restart);
        }
      else
        {
          local Vec u[m];       /* aliases to subsections */
          vec_aliases_create1 (U, m, u);

          if (bgy3d_getopt_test ("--load-guess"))
            bgy3d_vec_read1 ("u%d.bin", m, u);
          else
            {
              /* The very first iteration! */
              for (int i = 0; i < m; i++)
                {
                  const real eps = 80.0; /* FIXME: literal */
                  VecCopy (uc, u[i]);
                  VecScale (u[i], - ((eps - 1) / eps) * solvent[i].charge);
                }
            }
          vec_aliases_destroy1 (U, m, u);
        }

      {
        /*
          Work area for iterate_u().   Enclosing block is to make sure
          Ctx ctx is  only used by iterate_u() and  noone else. Though
          beware that the assignments  below just add other references
          to the local variables accessible by their name too.
        */
        Ctx ctx =
          {
            .BHD = BHD,                 /* in, mostly */
            .m = m,                     /* in */
            .solvent = solvent,         /* in */
            .rhos = rhos,               /* in */
            .kernel_fft = (void*) kernel_fft, /* in */
            .omega_fft = (void*) omega_fft,   /* in */
            .u0 = u0,                         /* in */
            .uc = uc,                         /* in */
            .x_lapl = x_lapl,                 /* inout */
            .g = g,                           /* out */
            .g_fft = g_fft,                   /* work */
            .du_acc_fft = du_acc_fft,         /* work */
            .work = work,                     /* work */
          };

        /*
          Find such a u that du as returned by iterate_u (&ctx, u, du)
          is zero. Cast  is there to silence the  mismatch in the type
          of  first   pointer  argument:  Ctx*  vs.    void*.   A  few
          alternative solvers  can be used here,  approximately in the
          order of sophistication (not preformace!):

          bgy3d_snes_picard(),
          bgy3d_snes_jager(),
          bgy3d_snes_newton().
        */
        bgy3d_snes_jager (PD, &ctx, (VectorFunc) iterate_u, U);
      }

      /*
        On  request,  save  u[]   to  binary  file.   If  you  provide
        --load-guess  in the  next  run this  u[]  can be  used as  an
        initial guess.  See above.
      */
      if (bgy3d_getopt_test ("--save-guess"))
        {
          local Vec u[m];       /* aliases to subsections */
          vec_aliases_create1 (U, m, u);

          bgy3d_vec_save1 ("u%d.bin", m, u);

          vec_aliases_destroy1 (U, m, u);
        }
    } /* for (damp = ... ) */


  /* Do not destroy  the long Vec U, return it as  restart info if the
     caller requested that: */
  if (restart)
    {
      /*
        This is probably QM code  calling, eventually this will not be
        the last call  during SCF.  Pass the PMF Vec  U (the long Vec)
        to help us restart iterations  in the next SCF round.  At this
        point we  give up ownersip  of the long  Vec U to  the caller.
        The caller will eventually have to dispose of it, unless it is
        passed back to this function.
      */
      *restart = (void*) U;
      U = NULL;                 /* because declared local */
    }
  else
    vec_pack_destroy1 (&U);     /* not vec_destroy()! */


  /*
    Compute, and eventually  (when Context** v is not  NULL) return to
    the  caller  the  iterator  over electrostatic  potential  of  the
    solvent. As  the distribution functions  are only used  to compute
    the  charge density  one can  supply either  g or  h =  g -  1 for
    neutral solvents.
  */
  Context *ret = info (BHD, m, solvent, n, solute, g, uc, uc_rho);
  if (medium)
    *medium = ret;
  else
    bgy3d_pot_destroy (ret);

  /* Clean up and exit ... */
  vec_destroy1 (m, u0);
  vec_destroy1 (m, g_fft);
  vec_destroy1 (m, x_lapl);

  vec_destroy (&du_acc_fft);
  vec_destroy (&work);
  vec_destroy (&uc);
  vec_destroy (&uc_rho);
}


/*
  Prepares the solvent (pair) properties  and calls the solver for the
  siglet distribution funcitons  g[].  Allocates g[m], deallocation is
  delegated to the caller.
*/
void bgy3d_solute_solve (const ProblemData *PD,
                         int m, const Site solvent[m],
                         int n, const Site solute[n],
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         Vec g[m],          /* out */
                         Context **medium,  /* out, optional */
                         Restart **restart) /* inout, optional */
{
  /* Show solvent/solute parameters: */
  bgy3d_sites_show ("Solvent", m, solvent);
  bgy3d_sites_show ("Solute", n, solute);

  PetscPrintf (PETSC_COMM_WORLD, "Solving BGY3dM (%d-site) equation ...\n", m);

  State *BHD = bgy3d_state_make (PD);

  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  /* Pair quantities  here, use symmetry wrt  (i <-> j)  to save space
     and work: */
  local Vec g2[m][m];               /* solvent-solvent pair distributions */
  vec_create2 (BHD->da, m, g2);

  /*
    Get g2[][] e.g. from a  previous pure solvent calculation. For CS2
    the original version hard-coded  reading the pair distributions as
    radial functions  from the  text files named  g2C, g2S,  and g2CS.
    This version uses g00.txt, g11.txt, and g01.txt instead.
  */
  if (bgy3d_getopt_test ("--from-radial-g2"))
    bgy3d_vec_read_radial2 (BHD->da, BHD->PD, "g%d%d.txt", m, g2);
  else
    bgy3d_vec_read2 ("g%d%d.bin", m, g2);

  local Vec omega_fft[m][m];    /* diagonal will be NULL */
  bgy3d_omega_fft_create (BHD, m, solvent, omega_fft); /* creates them */

  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)
      VecScale (omega_fft[i][j], PD->h[0] * PD->h[1] * PD->h[2]); /* historically */

  /*
    These  are the solvent  kernels, e.g.   HH, HO,  OH, OO  stored as
    complex vectors in momentum  space.  Note that kernel_fft[i][j] ==
    kernel_fft[j][i].
  */
  local Vec kernel_fft[m][m];
  vec_create2 (BHD->dc, m, kernel_fft); /* complex */

  /*
    Returns div (F * g2).  Note the calculation of F is divided due to
    long  range Coulomb  interation.   See comments  in the  function.
    Here F is force within solvents particles.

    The pairwise long-range interaction is included in the kernel.
  */
  solvent_kernel (BHD, m, solvent, g2, kernel_fft);

  /* Here the  storage for  the output is  allocated, the  caller will
     have to destroy them: */
  vec_create1 (BHD->da, m, g); /* real */

  solute_solve (BHD,
                m, solvent, kernel_fft, omega_fft, /* in */
                n, solute, density,                /* in */
                g, medium,                         /* out */
                restart);                          /* inout */

  /* Clean up and exit. Pair quantities here: */
  vec_destroy2 (m, g2);
  vec_destroy2 (m, kernel_fft);

  /* Diagonal is NULL: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)
      vec_destroy (&omega_fft[i][j]);

  /* Delegated to the caller: vec_destroy1 (m, g); */

  bgy3d_state_destroy (BHD);
}


/* This one emulates historical solver interface: */
Vec BGY3d_solute_solve (const ProblemData *PD, Vec g_ini)
{
  (void) g_ini;                 /* FIXME: interface obligation */

  int n;                        /* number of solute sites */
  const Site *solute;           /* solute[n] */

  int m;                        /* number of solvent sites */
  const Site *solvent;          /* solvent[m] */

  /* Get the number of solvent sites and their parameters: */
  bgy3d_solvent_get (&m, &solvent);

  char name[200] = "hydrogen chloride"; /* default solute */

  /* Solutes name, HCl by default: */
  bgy3d_getopt_string ("--solute", name, sizeof(name));

  /* Code used to be verbose: */
  PetscPrintf (PETSC_COMM_WORLD, "Solute is %s.\n", name);

  /* Get the solute from the tables: */
  bgy3d_solute_get (name, &n, &solute);

  /*
    This does  the real work. Vec  g[m] is intent(out)  in all senses,
    dont  forget   to  destroy   them.   Here  no   additional  charge
    distribution, so supply NULL for the function pointer:
  */
  Vec g[m];
  bgy3d_solute_solve (PD, m, solvent, n, solute,
                      NULL,     /* no electron density */
                      g,        /* out */
                      NULL,     /* dont need medium info */
                      NULL);    /* not going to restart */

  /* Save final distribution, use binary format: */
  bgy3d_vec_save1 ("g%d.bin", m, g);

  vec_destroy1 (m, g);

  return PETSC_NULL;            /* fake, interface obligation */
}

#ifdef WITH_EXTRA_SOLVERS
/*
  Deprecated!

  Convert   the   solvers  to   use   the   kernel   as  prepared   by
  solvent_kernel() and  applied by apply(), instead  of using weighted
  forces.

  Given  pairwise  distributions g2[][]  compute  weighted forces  and
  coulomb interction.

  Side  effects:  by  way  of  pair() uses  BHD->fg2_fft[3]  as  work
  arrays.
*/
void RecomputeInitialFFTs (State *BHD,
                           int m,
                           Vec g2[m][m],        /* real, in */
                           Vec fs_g2_fft[m][m][3], /* complex, out */
                           Vec fl_g2_fft[m][m][3], /* complex, out */
                           Vec u2[m][m],     /* real, out */
                           Vec u2_fft[m][m]) /* complex, out */
{
  assert (m == 2);              /* FIXME: uses global solvent[2] */

  PetscPrintf (PETSC_COMM_WORLD, "Recomputing FFT data\n");

  /* FIXME: maybe use BHD->v[3] for one of them? */
  Vec force_short[3];           /* work vectors for pair() */
  Vec force_long[3];            /* work vectors for pair() */
  vec_create1 (BHD->da, 3, force_short);
  vec_create1 (BHD->da, 3, force_long);


  /* Over all distinct solvent site pairs: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /* Does real work: */
        pair (BHD, solvent[i], solvent[j],
              g2[i][j],
              force_short, force_long, /* work vectors*/
              fs_g2_fft[i][j], fl_g2_fft[i][j],
              u2[i][j], u2_fft[i][j]); /* ij = ji */
      }

  /* Clean up and exit: */
  vec_destroy1 (3, force_short);
  vec_destroy1 (3, force_long);
}

/*
  Deprecated!

  Convert   the   solvers  to   use   the   kernel   as  prepared   by
  solvent_kernel() and  applied by apply(), instead  of using weighted
  forces.

  Side effects: none, but consider efficiency.
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
  Vec kernel_fft = vec_create (BHD->dc);
  Vec g_fft = vec_create (BHD->dc);
  Vec dg_fft = vec_create (BHD->dc);

  /*
    FIXME:  Move   computation  of  the  kernel  out   of  the  BGY3dM
    iterations.   Here we put  the fft  of the  kernel into  local Vec
    kernel_fft:
  */
  kernel (BHD->dc, BHD->PD, fg2_fft, coul_fft,
          kernel_fft);          /* result */

  /* fft(g) */
  MatMult (BHD->fft_mat, g, g_fft);

  /* Will be incremented: */
  VecSet (dg_fft, 0.0);

  /* Apply the kernel, Put result into the complex temp Vec dg_fft: */
  apply (kernel_fft, g_fft, scale, dg_fft);

  /* ifft(dg) */
  MatMultTranspose (BHD->fft_mat, dg_fft, dg);

  vec_destroy (&kernel_fft);
  vec_destroy (&g_fft);
  vec_destroy (&dg_fft);
}

/*
  Deprecated!

  Convert   the   solvers  to   use   the   kernel   as  prepared   by
  solvent_kernel() and  applied by apply(), instead  of using weighted
  forces.

  Side effects: none, but see Compute_H2O_interS_C().
 */
void Compute_H2O_interS (const State *BHD, /* NOTE: modifies BHD->fft dynamic arrays */
                         Vec fg2_fft[3],   /* complex, intent(in) */
                         Vec g,            /* real, intent(in) */
                         real rho, Vec dg_help)
{
    Compute_H2O_interS_C(BHD, fg2_fft, g, NULL, rho, dg_help);
}

/*
  In the  original code these  were the essential  differences between
  BGY3dM_solve_H2O_2site()  and BGY3dM_solve_H2O_3site().   First, the
  density of the H-sites was doubled:

  > BHD->rho_H = 2.0 * BHD->rho_H;

  Second,  there was an  additional "normalization"  term for  the H-H
  pair added to the δu accumulator for gH (called dg_new here):

  > Solve_NormalizationH2O_smallII (BHD, gH, r_HH, gH, tH , dg_new2, f, zpad);
  > Compute_dg_H2O_intra_ln (BHD, tH, r_HH, dg_new2, f);
  > VecAXPY (dg_new, 1.0, dg_new2);

  Third, the  corresponding O-H "normalization" term was  added to the
  δu accumulator for gO with a factor 2.0 as opposed to 1.0:

  < VecAXPY(dg_new, 1.0, dg_new2);
  ---
  > VecAXPY(dg_new, 2.0, dg_new2);

  Otherwise the code pf BGY3dM_solve_H2O_3site() appeared to be a copy
  of BGY3dM_solve_H2O_2site() with  a few inessential differences like
  using  a   different  default  solute  (butanoic   acid  for  3-site
  vs. hexane  for two-site) and  a different definition of  the mixing
  factor a0 = 0.1/(count+5.0) in the damp-loop.

  This approach appears to assume that the distributions of the H1 and
  H2  sites  around the  solute  is  the  same.  Probably  a  naturall
  assumption given  the symmetry. There is nothing  that would enforce
  this symmetry if treating the two sites idependently.
*/
Vec BGY3dM_solve_H2O_3site(const ProblemData *PD, Vec g_ini)
{
  (void) g_ini;                 /* FIXME: interface obligation */

  real a1, a, damp, damp_LJ;
  real count = 0.0;
  int iter;
  Vec dgH, dgO,  dg_new, dg_new2;
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

  /* Site specific  density.  Computed as a solvent  density rho times
     number of sites of that type in a solvent: */
  real rhos[m];
  for (int i = 0; i < m; i++)
    rhos[i] = PD->rho;
  rhos[0] = 2.0 * rhos[0];

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

  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem across iterations. See call to KSPSetInitialGuessNonzero()
    in bgy3d-poisson.c:
  */
  Vec x_lapl[2];                /* real */
  vec_create1 (BHD->da, 2, x_lapl);
  for (int i = 0; i < 2; i++)
    VecSet (x_lapl[i], 0.0);
#endif

  g[0] = vec_create (BHD.da);
  g[1] = vec_create (BHD.da);
  dgH = vec_create (BHD.da);
  dgO = vec_create (BHD.da);
  dg_new = vec_create (BHD.da);
  dg_new2 = vec_create (BHD.da);
  work = vec_create (BHD.da);

  tH = vec_create (BHD.da);
  tO = vec_create (BHD.da);

  dg_newH = vec_create (BHD.da);
  dg_newO = vec_create (BHD.da);

  uc = vec_create (BHD.da); /* common for all sites */

  dg_histO = vec_create (BHD.da);
  dg_histH = vec_create (BHD.da);
  VecSet(dg_histH, 0.0);
  VecSet(dg_histO, 0.0);

  Vec g0[2];
  for (int i = 0; i < 2; i++)
    g0[i] = vec_create (BHD.da);

  /* set initial guess*/
  VecSet(dgH,0);
  VecSet(dgO,0);
  VecSet(dg_new,0.0);

  /* load initial configuration from file ??? */
  if (bgy3d_getopt_test ("--load-guess")) {
      PetscPrintf(PETSC_COMM_WORLD,"Loading binary files...");
      dgH = bgy3d_vec_load ("dg0.bin"); /* dgH */
      dgO = bgy3d_vec_load ("dg1.bin"); /* dgO */
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

      RecomputeInitialFFTs (&BHD, 2,
                            g2,            /* real, in */
                            BHD.fs_g2_fft, /* complex, out */
                            BHD.fl_g2_fft, /* complex, out */
                            BHD.u2,        /* real, out */
                            BHD.u2_fft);   /* complex, out */


      /* Compute solute field in this block:*/
      {
        int n;                        /* number of solute sites */
        const Site *sites;            /* [n], array of sites */

        /* Get the solute from the tables: */
        bgy3d_solute_get ("butanoic acid", &n, &sites);

        /* This does the real work: */
        bgy3d_solute_field (&BHD,
                            2, solvent,
                            g0, uc, /* intent(out) */
                            n, sites,
                            NULL); /* void (*density)(...) */
      }

      PetscPrintf(PETSC_COMM_WORLD,"New lambda= %f\n", a0);

      /* Historically   short-range  potential   is   stored  with   a
         factor: */
      for (int i = 0; i < 2; i++)
        VecScale (g0[i], BHD.PD->beta);

      bgy3d_impose_laplace_boundary (&BHD, g0[0], x_lapl[0]);
      bgy3d_impose_laplace_boundary (&BHD, g0[1], x_lapl[1]);

      /* g=g0*exp(-dg) */
      bgy3d_compute_g (g[0], g0[0], dgH);
      bgy3d_compute_g (g[1], g0[1], dgO);

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
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[0][1], g[1], rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[0][0], g[0], rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][1], g[1], BHD.u2_fft[0][1], damp * rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][0], g[0], BHD.u2_fft[0][0], damp * rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);

          /* Specific to H2O: */
          Solve_NormalizationH2O_smallII (&BHD, g[0], r_HH, g[0], tH , dg_new2, work);

          Compute_dg_H2O_intra_ln(&BHD, tH, r_HH, dg_new2);
          VecCopy (dg_new2, work); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);
          /* end specific */

          Solve_NormalizationH2O_smallII (&BHD, g[0], r_HO, g[1], tO , dg_new2, work);

          Compute_dg_H2O_intra_ln(&BHD, tO, r_HO, dg_new2);
          VecCopy (dg_new2, work); /* FIXME: need that? */
          VecAXPY(dg_new, 1.0, dg_new2);

          /* Copy  the electrostatic potential  scaled by  the solvent
             site charges into predefined locations: */
          VecAXPY(dg_new, solvent[0].charge, uc);

          bgy3d_impose_laplace_boundary (&BHD, dg_new, x_lapl[0]);

          VecCopy(dg_new, dg_newH);

          /* O */
          VecSet(dg_new,0.0);
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[1][1], g[1], rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS(&BHD, BHD.fs_g2_fft[0][1], g[0], rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          VecScale(dg_new,damp_LJ);

          /* Coulomb long */
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[1][1], g[1], BHD.u2_fft[1][1], damp * rhos[1], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);
          Compute_H2O_interS_C(&BHD, BHD.fl_g2_fft[0][1], g[0], BHD.u2_fft[0][1], damp * rhos[0], dg_new2);
          VecAXPY(dg_new, 1.0, dg_new2);

          Solve_NormalizationH2O_smallII (&BHD, g[1], r_HO, g[0], tH , dg_new2, work);
          Compute_dg_H2O_intra_ln(&BHD, tH, r_HO, dg_new2);
          VecCopy (dg_new2, work); /* FIXME: need that? */

          VecAXPY(dg_new, 2.0, dg_new2); /* factor 2.0 specific to H2O */

          /* Copy  the electrostatic potential  scaled by  the solvent
             site charges into predefined locations: */
          VecAXPY(dg_new, solvent[1].charge, uc);

          bgy3d_impose_laplace_boundary (&BHD, dg_new, x_lapl[1]);

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
          bgy3d_compute_g (g[0], g0[0], dgH);
          bgy3d_compute_g (g[1], g0[1], dgO);

          PetscPrintf (PETSC_COMM_WORLD, "Q=% e ", ComputeCharge (PD, 2, solvent, g));

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
      bgy3d_vec_save_ascii (nameH, g[0]); /* g_H */
      bgy3d_vec_save_ascii (nameO, g[1]); /* g_O */
      PetscPrintf(PETSC_COMM_WORLD,"done.\n");

      /* save g to binary file */
      if (bgy3d_getopt_test ("--save-guess")) {
          PetscPrintf(PETSC_COMM_WORLD,"Writing binary files...");
          bgy3d_vec_save ("dg0.bin", dgH); /* dgH */
          bgy3d_vec_save ("dg1.bin", dgO); /* dgO */
          PetscPrintf(PETSC_COMM_WORLD,"done.\n");
      }
    }

  vec_destroy1 (2, g);
  vec_destroy1 (2, g0);

  vec_destroy (&dgH);
  vec_destroy (&dgO);
  vec_destroy (&dg_new);
  vec_destroy (&dg_new2);
  vec_destroy (&work);

  vec_destroy (&tH);
  vec_destroy (&tO);

  vec_destroy (&uc);

  vec_destroy (&dg_newH);
  vec_destroy (&dg_newO);
  vec_destroy (&dg_histH);
  vec_destroy (&dg_histO);

  finalize_state (&BHD);

  return PETSC_NULL;
}
#endif
