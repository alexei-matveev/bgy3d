/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2OS.c,v 1.20 2007-08-03 15:59:50 jager Exp $ */
/*==========================================================*/

/*
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
  the solution of u can be repsented by a difference of two functions:

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

  and sum  to the solution  by the end  of each iteration  for solvent
  site A and B, in the code hereafter.  These site specific potentials
  can  be  obtained  by  multiplying  the  common  long-range  Coulomb
  potential (stored in Vec uc) by the solvent site charge.
*/


#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* needs Site */
#include "bgy3d-force.h"        /* Coulomb_short_grad() */
#include "bgy3d-getopt.h"
#include "bgy3d-vec.h"
#include "bgy3d-pure.h"
#include "bgy3d-poisson.h"      /* bgy3d_poisson() */
#include "bgy3d-dirichlet.h"    /* Laplace staff */
#include "bgy3d-potential.h"    /* Context, etc. */
#include "bgy3d-fftw.h"         /* bgy3d_fft_interp() */
#include "bgy3d-snes.h"         /* bgy3d_snes_newton() */
#include "bgy3d-impure.h"

#ifndef L_BOUNDARY_MG
#include "bgy3d-multigrid.h"    /* InitializeDMMGSolver */
#endif

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
  bgy3d_pair (BHD, a, b, f_short, f_long, NULL, NULL, u2, u2_fft,
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
  Vec work = BHD->v[0];

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
  int x[3], n[3], i[3], N[3];

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
  DAVecRestoreArray (dc, dfg, &dfg_);
  if (coul)
    DAVecRestoreArray (dc, coul, &coul_);
  FOR_DIM
    DAVecRestoreArray (dc, fg[dim], &fg_[dim]);

  /*
    Translate  the kernel  so  that the  real-space representation  is
    centered at the  grid corner (0, 0, 0) and not  at the grid center
    like other  grid representations. This  is what is assumed  in the
    convolution integrals:
  */
  bgy3d_vec_fft_trans (dc, N, dfg);
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

  /*
    Retrive the  precomuted Fourier transform of of  the divergence of
    the  "weighted  force"  vector  div   (F  g)  which  serves  as  a
    convolution  kernel in  BGY3dM equations.   The  additional factor
    -k^(-2) effectively  included in  the kernel recovers  the Fourier
    transform of the corresponding Poisson solution.  The factor scale
    = βρ is not included, on the other hand. See kernel() for details:

    Loop over local portion of grid:
  */
  for (int k = x[2]; k < x[2] + n[2]; k++)
    for (int j = x[1]; j < x[1] + n[1]; j++)
      for (int i = x[0]; i < x[0] + n[0]; i++)
        dg_[k][j][i] += scale * (ker_[k][j][i] * g_[k][j][i]); /* complex */

  DAVecRestoreArray (dc, ker, &ker_);
  DAVecRestoreArray (dc, dg, &dg_);
  DAVecRestoreArray (dc, g, &g_);
}

/*
  Given  the  solvent  description  and  pair  distribution  functions
  compute the multiplicative  complex momentum-space representation of
  kernel to be applied every iteration.

  Side  effects:  by  way  of  pair() uses  BHD->fg2_fft[3]  as  work
  arrays.
*/
static void solvent_kernel (State *BHD, int m, const Site solvent[m],
                            Vec g2[m][m],      /* real, in */
                            Vec ker_fft[m][m]) /* complex, out */
{
  /* Real  work  vectors,  re-used   for  all  m(m+1)/2  solvent  site
     pairs. */
  Vec u2 = bgy3d_vec_create (BHD->da);
  Vec fs[3], fl[3];
  bgy3d_vec_create1 (BHD->da, 3, fs); /* 3-vector */
  bgy3d_vec_create1 (BHD->da, 3, fl); /* 3-vector */

  /* Complex work vectors, re-used for all pairs: */
  Vec kl_fft = bgy3d_vec_create (BHD->dc);
  Vec u2_fft = bgy3d_vec_create (BHD->dc);
  Vec fs_g2_fft[3], fl_g2_fft[3];
  bgy3d_vec_create1 (BHD->dc, 3, fs_g2_fft); /* 3-vector */
  bgy3d_vec_create1 (BHD->dc, 3, fl_g2_fft); /* 3-vector */

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
        kernel (BHD->dc, BHD->PD, fs_g2_fft, NULL, ker_fft[i][j]);

        /* This used to be  the long-range contribution, we compute it
           separately, but later add to the total kernel: */
        kernel (BHD->dc, BHD->PD, fl_g2_fft, u2_fft, kl_fft);

        /*
          FIXME:  what is  the point  to split  the kernel  in two
          pieces?  Redefine S := S + L and forget about L:

          ker = ker_l + ker_s
        */
        VecAXPY (ker_fft[i][j], 1.0, kl_fft);
      }

  /* Clean up and exit: */
  bgy3d_vec_destroy (&u2);
  bgy3d_vec_destroy (&u2_fft);
  bgy3d_vec_destroy (&kl_fft);
  bgy3d_vec_destroy1 (3, fs);
  bgy3d_vec_destroy1 (3, fl);
  bgy3d_vec_destroy1 (3, fs_g2_fft);
  bgy3d_vec_destroy1 (3, fl_g2_fft);
}

/* Dipole of the cores of a single specie: */
static void dipole (int n, const Site sites[n], real d[3], real *d_norm)
{
  FOR_DIM
    d[dim] = 0.0;

  for (int i = 0; i < n; i++)
    FOR_DIM
      d[dim] += sites[i].charge * sites[i].x[dim];

  *d_norm = sqrt (SQR (d[0]) + SQR (d[1]) + SQR (d[2]));
}

static void moments (const State *BHD, Vec v,
                     real *q, real *x, real *y, real *z)
{
  const real *h = BHD->PD->h;   /* h[3] */
  const real h3 = h[0] * h[1] * h[2];

  bgy3d_vec_moments (BHD->da, v, q, x, y, z);

  *q *= h3 * 1;
  *x *= h3 * h[0];
  *y *= h3 * h[1];
  *z *= h3 * h[2];
}

static real ComputeCharge (const ProblemData *PD,
                           int m, const Site solvent[m],
                           const Vec g[m])
{
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];

  real total = 0.0;
  for (int i = 0; i < m; i++)
    total += solvent[i].charge * bgy3d_vec_sum (g[i]);

  /* FIXME: different rhos[]? */
  return total * h3 * PD->rho;
}

/*
  Function to  generate solvent field that  is supposed to  act on the
  solute  electrons. At  the moment  only the  electrostatic  field is
  computed.   This  "reaction" field  should  be  consistent with  the
  "action" of  the solute electrons  on the solvent sites  computed in
  bgy3d-solvent.c  that  is also  pure  electrostatics. Returns  both,
  electrostatic potential  with a  boundary correction and  the actual
  solvent charge density.
*/
static void bgy3d_solvent_field (const State *BHD, /* intent(in) */
                                 int m, const Site solvent[m],
                                 Vec g[m],        /* intent(in) */
                                 Vec ve, Vec rho) /* intent(out) */
{
  /* FIXME: this  code assumes  the same density  rho for  all solvent
     particles. */

  /*
    Solvent   charge  density   for  Poisson   equation   and  (later)
    integration   with   solute   electrostatic  potential   will   be
    accumulated here:
  */
  VecSet (rho, 0.0);
  for (int i = 0; i < m; i++)
    VecAXPY (rho, solvent[i].charge * BHD->PD->rho, g[i]);

  /*
    Solve Poisson equation for rho.  Note that the output potential is
    in kcals  (see the  definiton of EPSILON0INV)  as all  energies in
    this code are:
  */
  bgy3d_poisson (BHD, ve, rho, -4 * M_PI * EPSILON0INV);

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
    Vec x = bgy3d_vec_create (BHD->da);

    bgy3d_impose_laplace_boundary (BHD, ve, x);

    bgy3d_vec_destroy (&x);
  }

  bgy3d_vec_save ("ve.bin", ve); /* for debugging only */
}

static void print_table (int n, const Site sites[n], const real vs[n])
{
  /* Average  over   all  sites,  has   no  real  meaning.   Only  for
     presentation purposes: */
  real v_avg = 0.0;
  for (int i = 0; i < n; i++)
    v_avg += vs[i];
  v_avg /= n;

  PetscPrintf (PETSC_COMM_WORLD,
               "#\t site\t x        \t y        \t z        \t q        \t δv        \t v0\n");
  for (int i = 0; i < n; i ++)
    PetscPrintf (PETSC_COMM_WORLD,
                 "%d\t%5s\t% f\t% f\t% f\t% f\t% f\t% f\n",
                 i + 1, sites[i].name,
                 sites[i].x[0], sites[i].x[1],  sites[i].x[2],
                 sites[i].charge, vs[i] - v_avg, v_avg);
}


/* g := exp (-u) */
static void mexp (Vec g, Vec u)
{
  real pure f (real u)
  {
    return exp (-u);
  }
  bgy3d_vec_map1 (g, f, u);
}


static void iterate (State *BHD,
                     int m,
                     const Site solvent[m],
                     const real rhos[m],
                     Vec ker_fft[m][m], /* in */
                     Vec omega[m][m],   /* in */
                     Vec u0[m],         /* in */
                     Vec uc,            /* in */
                     Vec u[m],          /* in */
                     Vec x_lapl[m],     /* inout */
                     Vec g[m],          /* out */
                     Vec g_fft[m],      /* work */
                     Vec du_acc_fft,    /* work */
                     Vec work,          /* work */
                     Vec du[m])         /* out */
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
        apply (BHD->dc, ker_fft[i][j], g_fft[j], beta * rhos[j],
               du_acc_fft); /* incremented! */

      /*
        Compute  IFFT  of  du_acc_fft  for  the  current  site.  Other
        contributions are added to the real space du_acc below:
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
          bgy3d_nssa_intra_log (BHD, g_fft[i], omega[i][j], g[j], work);

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
  Vec *ker_fft;                 /* [m][m] */
  Vec *omega;                   /* [m][m] */
  Vec *u0;                      /* [m] */
  Vec uc;                       /*  */
  Vec *x_lapl;                  /* [m] */
  Vec *g;                       /* [m] */
  Vec *g_fft;                   /* [m] */
  Vec du_acc_fft;               /* work */
  Vec work;                     /* work */
} Ctx;


static void iterate_u (Ctx *s, Vec us, Vec dus)
{
  const int m = s->m;           /* number of solvent sites */

  /* Compiler expects us to pass [m][m] arrays to iterate(): */
  Vec (*ker_fft)[m] = (void*) s->ker_fft;
  Vec (*omega)[m] = (void*) s->omega;

  /* Establish  aliases to  the subsections  of  the long  Vec us  and
     dus */
  Vec u[m], du[m];
  bgy3d_vec_aliases_create (us, m, u);
  bgy3d_vec_aliases_create (dus, m, du);

  /* Iterate u[] -> du[]: */
  iterate (s->BHD,              /*  1 in */
           s->m,                /*  2 in */
           s->solvent,          /*  3 in */
           s->rhos,             /*  4 in */
           ker_fft,             /*  5 in */
           omega,               /*  6 in */
           s->u0,               /*  7 in */
           s->uc,               /*  8 in */
           u,                   /*  9 in <- here */
           s->x_lapl,           /* 10 inout */
           s->g,                /* 11 out */
           s->g_fft,            /* 12 work */
           s->du_acc_fft,       /* 13 work */
           s->work,             /* 14 work */
           du);                 /* 15 out <- here */

  /* This  destroys the  aliases, but  does  not free  the memory,  of
     course. The actuall data is owned by Vec us and Vec dus: */
  bgy3d_vec_aliases_destroy (m, u);
  bgy3d_vec_aliases_destroy (m, du);
}

/*
  This function is the main entry  point for the BGY3dM equation for a
  m-site solvent and an arbitrary solute.  The vectors in

  Vec g[m], intent(out)

  are  initialized as global  distributed arrays  and filled  with the
  solvent site  distributions. It is the responsibility  of the caller
  to destroy them when no more needed.

  Context **v, inent(out);

  Is set to  an iterator over the solvent  potential field. The caller
  is responsible  for calling bgy3d_pot_destroy() when it  is not more
  needed.
*/
void bgy3d_solute_solve (const ProblemData *PD,
                         int m, const Site solvent[m],
                         int n, const Site solute[n],
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         Vec g[m],    /* intent(out) */
                         Context **v) /* intent(out) */
{
  /* Show solvent/solute parameters: */
  bgy3d_sites_show ("Solvent", m, solvent);
  bgy3d_sites_show ("Solute", n, solute);

  int namecount = 0;

  PetscPrintf (PETSC_COMM_WORLD, "Solving BGY3dM (%d-site) equation ...\n", m);

  State *BHD = bgy3d_state_make (PD);

  /* Code used to be verbose: */
  bgy3d_state_print (BHD);

  if (r_HH > 0.0)
    PetscPrintf (PETSC_COMM_WORLD, "WARNING: Solvent not a 2-Site model!\n");

  /* Site specific  density.  Computed as a solvent  density rho times
     number of sites of that type in a solvent: */
  real rhos[m];
  for (int i = 0; i < m; i++)
    rhos[i] = PD->rho;          /* 2 * PD->rho; */

  /* Here dN is rho * dV, with dV being the weight of a grid point: */
  const real dN = PD->rho * PD->h[0] * PD->h[1] * PD->h[2];
  {
    const int N3 =  PD->N[0] * PD->N[1] * PD->N[2];
    PetscPrintf (PETSC_COMM_WORLD, "Number of solvent molecules is %f\n", N3 * dN);
  }

  /* Pair quantities  here, use symmetry wrt  (i <-> j)  to save space
     and work: */
  Vec g2[m][m];               /* solvent-solvent pair distributions */
  bgy3d_vec_create2 (BHD->da, m, g2);

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

  /*
    These  are the solvent  kernels, e.g.   HH, HO,  OH, OO  stored as
    complex  vectors in  momentum space.   Note that  ker_fft[i][j] ==
    ker_fft[j][i].
  */
  Vec ker_fft[m][m];
  bgy3d_vec_create2 (BHD->dc, m, ker_fft); /* complex */

  /*
    Returns div (F * g2).  Note the calculation of F is divided due to
    long  range Coulomb  interation.   See comments  in the  function.
    Here F is force within solvents particles.

    The pairwise long-range interaction is included in the kernel.
  */
  solvent_kernel (BHD, m, solvent, g2, ker_fft);

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
  Vec x_lapl[m];                /* real */
  bgy3d_vec_create1 (BHD->da, m, x_lapl);
  for (int i = 0; i < m; i++)
    VecSet (x_lapl[i], 0.0);
#endif

#ifdef L_BOUNDARY_MG
  /* Create KSP environment */
  InitializeDMMGSolver(BHD);
#endif

  Vec omega[m][m];
  {
    /* FIXME: m  x m  distance matrix does  not handle  equivalent sites
       well.  Diagonal zeros are never referenced: */
    real r[m][m];
    bgy3d_sites_dist_mat (m, solvent, r);

    /* Precompute omega[][]: */
    for (int i = 0; i < m; i++)
      for (int j = 0; j < i; j++)
        {
          omega[j][i] = omega[i][j] = bgy3d_vec_create (BHD->dc);
          bgy3d_omega (BHD->PD, BHD->dc, r[i][j], omega[i][j]);
        }
  }

  /*
    These complex  vectors will  hold FFT of  the current  g. Allocate
    enough to  hold a  local portion  of the grid  and free  after the
    loop.
  */
  Vec g_fft[m];
  bgy3d_vec_create1 (BHD->dc, m, g_fft); /* complex */

  /* Here the  storage for  the output is  allocated, the  caller will
     have to destroy them: */
  bgy3d_vec_create1 (BHD->da, m, g); /* real */

  /*
    There is no  point to transform each contribution  computed in the
    momentum space back, accumulate them on the k-grid in this complex
    Vec:
  */
  Vec du_acc_fft = bgy3d_vec_create (BHD->dc); /* complex */

  Vec work = bgy3d_vec_create (BHD->da);

  /* Coulomb long, common for all sites: */
  Vec uc = bgy3d_vec_create (BHD->da);

  /* Electron density for integration: */
  Vec uc_rho = bgy3d_vec_create (BHD->da);

  /*
    Later u0  = beta *  (VM_LJ + VM_coulomb_short), which  is -log(g0)
    actually.  See: (5.106)  and (5.108) in Jager's thesis.  It is not
    filled with data yet, I assume.
  */
  Vec u0[m];                          /* real */
  bgy3d_vec_create1 (BHD->da, m, u0); /* real */

  /*
    For primary  variable there are two  ways to access  the data: via
    the long Vec and m shorter  Vecs aliased to the subsections of the
    longer one.
  */
  Vec us = bgy3d_vec_pack_create (BHD->da, m);  /* long Vec */
  Vec dus = bgy3d_vec_pack_create (BHD->da, m); /* long Vec */

  Vec u[m], du[m];         /* in- and output of the BGY3D iteration */
  bgy3d_vec_aliases_create (us, m, u);   /* aliases to subsections */
  bgy3d_vec_aliases_create (dus, m, du); /* aliases to subsections */

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
        Set initial guess, either here or by reading from file. At the
        end of the "dump" loop du[] is written to disk, so that in the
        next iteration we will read an updated version:
      */
      if (bgy3d_getopt_test ("--load-guess"))
        bgy3d_vec_read1 ("u%d.bin", m, u);
      else
        for (int i = 0; i < m; i++)
          VecSet (u[i], 0.0);

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
            .ker_fft = (void*) ker_fft, /* in */
            .omega = (void*) omega,     /* in */
            .u0 = u0,                   /* in */
            .uc = uc,                   /* in */
            .x_lapl = x_lapl,           /* inout */
            .g = g,                     /* out */
            .g_fft = g_fft,             /* work */
            .du_acc_fft = du_acc_fft,   /* work */
            .work = work,               /* work */
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
        bgy3d_snes_jager (PD, &ctx, (Function) iterate_u, us);
      }

      /* FIXME:  Debug  output  from  every iteration  with  different
         overall scale factors damp.  Remove when no more needed. */
      {
        char fmt[20];
        snprintf (fmt, sizeof fmt, "vec%%d-%d.m", namecount++);
        bgy3d_vec_save_ascii1 (fmt, m, g);
      }

      /* Save du to binary file. FIXME: Why du and not g? */
      if (bgy3d_getopt_test ("--save-guess"))
        bgy3d_vec_save1 ("u%d.bin", m, u);

    } /* for (damp = ... ) */

  /*
    Compute, and eveutually  (when Context** v is not  NULL) return to
    the  caller  the  iterator  over electrostatic  potential  of  the
    solvent:
  */
  {
    /* Solvent electrostaic potential field: */
    Vec ve = bgy3d_vec_create (BHD->da);

    /* Keep solvent charge density for integration: */
    Vec ve_rho = bgy3d_vec_create (BHD->da);

    /* This fills Vec ve with solvent electrostatic potential: */
    bgy3d_solvent_field (BHD, m, solvent, g, ve, ve_rho);

    /* Optionally, return the iterator over the solvent field: */
    if (v)
      *v = bgy3d_pot_create (BHD, ve);

    /* Also compute dipole moment: */
    real d[3], d_norm;
    dipole (n, solute, d, &d_norm);

    PetscPrintf (PETSC_COMM_WORLD,
                 "|<x|ρ_n>| = |% f, % f, % f| = %f (dipole moment of solute cores)\n",
                 d[0], d[1], d[2], d_norm);

    {
      real q, x, y, z;
      moments (BHD, ve_rho, &q, &x, &y, &z);

      const real d = sqrt (SQR (x) + SQR (y) + SQR (z));

      PetscPrintf (PETSC_COMM_WORLD,
                   "|<x|ρ_v>| = |% f, % f, % f| = %f (dipole moment of solvent medium)\n",
                   x, y, z, d);
      PetscPrintf (PETSC_COMM_WORLD, "q_v = %f (charge of solvent medium)\n", q);
    }
    /*
      Integration of:

      1.  Interaction  energy of solute  point cores with  the solvent
      electrostatic field.

      2. Solvent electrostatic field with diffuse solute electron
      density.

      3. Solvent charge density with long-range solute electrostatic
      field.

      Vec uc and uc_rho  are returned by bgy3d_solute_field().  Vec ve
      and ve_rho are obtained from bgy3d_solvent_field().
    */

    /* 1. */
    if (v)                 /* If it was requested by the caller ... */
      {
        real xs[n][3], vs[n];   /* coordinates and potential values */

        /* Extract coordinates into plain array: */
        for (int i = 0; i < n; i++)
          FOR_DIM
            xs[i][dim] = solute[i].x[dim];

        bgy3d_pot_interp (*v, n, xs, vs);

        print_table (n, solute, vs);

        real val1 = 0.0;
        for (int i = 0; i < n; i++)
          val1 += solute[i].charge * vs[i];

        PetscPrintf (PETSC_COMM_WORLD,
                     "<U_v|ρ_N> = %lf (solvent electrostatic field with solute point nuclei)\n",
                     val1);
      }

    /* 2 and 3. */
    real val2, val3;
    {
      /* Dot-product is an integral over space (up to a factor): */
      const real h3 = BHD->PD->h[0] * BHD->PD->h[1] * BHD->PD->h[2];
      val2 = h3 * bgy3d_vec_dot (ve, uc_rho);
      val3 = h3 * bgy3d_vec_dot (uc, ve_rho);
    }
    bgy3d_vec_destroy (&ve);            /* yes, we do! */
    bgy3d_vec_destroy (&ve_rho);

    PetscPrintf (PETSC_COMM_WORLD,
                 "<U_v|ρ_u> = %lf "
                 "(solvent electrostatic field with diffuse charge density of solute)\n",
                 val2);
    PetscPrintf (PETSC_COMM_WORLD,
                 "<ρ_v|U_u> = %lf "
                 "(solvent charge density with long-range electrostatic field of solute)\n",
                 val3);
    PetscPrintf (PETSC_COMM_WORLD,
                 "     diff = % lf\n", val3 - val2);
  }

  /* Clean up and exit ... */
  bgy3d_vec_aliases_destroy (m, u);
  bgy3d_vec_aliases_destroy (m, du);

  bgy3d_vec_pack_destroy (&us);  /* not bgy3d_vec_destroy()! */
  bgy3d_vec_pack_destroy (&dus); /* not bgy3d_vec_destroy()! */

  bgy3d_vec_destroy1 (m, u0);
  bgy3d_vec_destroy1 (m, g_fft);
  bgy3d_vec_destroy1 (m, x_lapl);

  /* Pair quantities here: */
  bgy3d_vec_destroy2 (m, g2);
  bgy3d_vec_destroy2 (m, ker_fft);

  bgy3d_vec_destroy (&du_acc_fft);
  bgy3d_vec_destroy (&work);
  bgy3d_vec_destroy (&uc);
  bgy3d_vec_destroy (&uc_rho);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)
      bgy3d_vec_destroy (&omega[i][j]);

  /* Delegated to the caller: bgy3d_vec_destroy1 (m, g); */

  bgy3d_state_destroy (BHD);
}

/* This one emulates historical solver interface: */
Vec BGY3dM_solve_H2O_2site (const ProblemData *PD, Vec g_ini)
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
  bgy3d_solute_solve (PD, m, solvent, n, solute, NULL, g, NULL);

  /* Save final distribution, use binary format: */
  bgy3d_vec_save1 ("g%d.bin", m, g);

  bgy3d_vec_destroy1 (m, g);

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
  bgy3d_vec_create1 (BHD->da, 3, force_short);
  bgy3d_vec_create1 (BHD->da, 3, force_long);


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
  bgy3d_vec_destroy1 (3, force_short);
  bgy3d_vec_destroy1 (3, force_long);
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
  Vec ker_fft = bgy3d_vec_create (BHD->dc);
  Vec g_fft = bgy3d_vec_create (BHD->dc);
  Vec dg_fft = bgy3d_vec_create (BHD->dc);

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

  bgy3d_vec_destroy (&ker_fft);
  bgy3d_vec_destroy (&g_fft);
  bgy3d_vec_destroy (&dg_fft);
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
  bgy3d_vec_create1 (BHD->da, 2, x_lapl);
  for (int i = 0; i < 2; i++)
    VecSet (x_lapl[i], 0.0);
#endif

#ifdef L_BOUNDARY_MG
  /* Create KSP environment */
  InitializeDMMGSolver(&BHD);
#endif

  g[0] = bgy3d_vec_create (BHD.da);
  g[1] = bgy3d_vec_create (BHD.da);
  dgH = bgy3d_vec_create (BHD.da);
  dgO = bgy3d_vec_create (BHD.da);
  dg_new = bgy3d_vec_create (BHD.da);
  dg_new2 = bgy3d_vec_create (BHD.da);
  work = bgy3d_vec_create (BHD.da);

  tH = bgy3d_vec_create (BHD.da);
  tO = bgy3d_vec_create (BHD.da);

  dg_newH = bgy3d_vec_create (BHD.da);
  dg_newO = bgy3d_vec_create (BHD.da);

  uc = bgy3d_vec_create (BHD.da); /* common for all sites */

  dg_histO = bgy3d_vec_create (BHD.da);
  dg_histH = bgy3d_vec_create (BHD.da);
  VecSet(dg_histH, 0.0);
  VecSet(dg_histO, 0.0);

  Vec g0[2];
  for (int i = 0; i < 2; i++)
    g0[i] = bgy3d_vec_create (BHD.da);

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

  bgy3d_vec_destroy1 (2, g);
  bgy3d_vec_destroy1 (2, g0);

  bgy3d_vec_destroy (&dgH);
  bgy3d_vec_destroy (&dgO);
  bgy3d_vec_destroy (&dg_new);
  bgy3d_vec_destroy (&dg_new2);
  bgy3d_vec_destroy (&work);

  bgy3d_vec_destroy (&tH);
  bgy3d_vec_destroy (&tO);

  bgy3d_vec_destroy (&uc);

  bgy3d_vec_destroy (&dg_newH);
  bgy3d_vec_destroy (&dg_newO);
  bgy3d_vec_destroy (&dg_histH);
  bgy3d_vec_destroy (&dg_histO);

  finalize_state (&BHD);

  return PETSC_NULL;
}
#endif
