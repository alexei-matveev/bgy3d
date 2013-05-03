/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3dH2O.c,v 1.42 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* struct Site */
#include "bgy3d-solvents.h"     /* needs Site */
#include "bgy3d-force.h"        /* bgy3d_force() */
#include "bgy3d-getopt.h"
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-dirichlet.h"    /* Laplace staff */
#include "bgy3d-vec.h"          /* vec_map*() */
#include "bgy3d-pure.h"
#include <complex.h>            /* after fftw.h */

static const real NORM_REG = 1.0e-1;
static const real NORM_REG2 = 1.0e-2;

/* FIXME: bgy3d-solvents.h pollutes the namespace: */
#undef sH
#undef eH
#undef qH
#undef sO
#undef eO
#undef qO


/* g := exp[-(u0 + du)], with a sanity check: */
static void compute_g (Vec g, Vec u0, Vec du)
{
  real pure f (real x, real y)
  {
    const real z = exp (- (x + y));
    assert (!isinf (z));
    assert (!isnan (z));
    return z;
  }
  vec_map2 (g, f, u0, du);
}


static void RecomputeInitialData (State *BHD,
                                  int m, const Site solvent[m],
                                  Vec f[m][m][3], Vec f_l[m][m][3],
                                  Vec u_ini[m][m], Vec c2[m][m],
                                  Vec u2[m][m], Vec u2_fft[m][m],
                                  real damp, real damp_LJ)
{
  PetscPrintf (PETSC_COMM_WORLD,
               "Recomputing initial data with damping factor %f (damp_LJ=%f)\n",
               damp, damp_LJ);

  /* Over all pairs: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        /*
          This  computes short-  and  long-rage forces,  corresponding
          potentials and more. Note  that we put the short-range force
          into BHD->F temporarily:
        */
        bgy3d_force (BHD, solvent[i], solvent[j],
                     f[i][j], f_l[i][j],
                     u_ini[i][j], c2[i][j],
                     u2[i][j], u2_fft[i][j],
                     damp, damp_LJ);

        /* Was previousely done in pair(): */
        VecScale (u_ini[i][j], BHD->PD->beta);

        /* BHD->F[][]  is a  total force,  apparently.  Long  range is
           stored also separately: */
        FOR_DIM
          VecAXPY (f[i][j][dim], 1.0, f_l[i][j][dim]);

      }
}

/* All vectors are complex here. No side effects. */
static void kapply (const State *BHD,
                    Vec fg2_fft[3],          /* intent(in) */
                    Vec g_fft, Vec coul_fft, /* intent(in) */
                    Vec dg_fft)              /* intent(out) */
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */
  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];
  const real fac = L / (2. * M_PI); /* BHD->f ist nur grad U, nicht F=-grad U  */

  complex ***g_fft_, ***dg_fft_, ***coul_fft_, ***fg2_fft_[3];
  DMDAVecGetArray (BHD->dc, g_fft, &g_fft_);
  DMDAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  DMDAVecGetArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DMDAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* Get local portion of the k-grid */
  int x[3], n[3], i[3];
  DMDAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);

          complex div = 0.0;    /* complex */
          if (likely (k2 != 0))
            {
              FOR_DIM
                div += ic[dim] * fg2_fft_[dim][i[2]][i[1]][i[0]];

              const real k_fac = (h3 * fac) / k2;

              /* "I" is an imaginary unit here: */
              div *= -I * k_fac;

              /* Long  range  Coulomb part.  Note  there  is not  "-I"
                 factor here: */
              div += coul_fft_[i[2]][i[1]][i[0]];

              div *=  h3 * g_fft_[i[2]][i[1]][i[0]];
            }
          dg_fft_[i[2]][i[1]][i[0]] = div;
        }
  DMDAVecRestoreArray (BHD->dc, g_fft, &g_fft_);
  DMDAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  DMDAVecRestoreArray (BHD->dc, coul_fft, &coul_fft_);
  FOR_DIM
    DMDAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /*
    Translate  the dg_fft  so  that the  real-space representation  is
    localized at  the grid center like other  grid representations and
    not at the grid corner.   This amounts to scaling the k-components
    by a phase factor:
  */
  bgy3d_vec_fft_trans (BHD->dc, N, dg_fft);
}

/* Side effects: uses one real and four complex work vectors.  Does (4
   + 1) FFTs. One inverse. */
static void Compute_dg_inter (State *BHD,
                              Vec fab_s[3], Vec fab_l[3], Vec gab,
                              Vec gb,
                              Vec cab_fft, /* long range coulomb */
                              real rhob,
                              Vec dua_fft,
                              Vec dua) /* intent(out) */
{
  const ProblemData *PD = BHD->PD;
  const real L = PD->interval[1] - PD->interval[0];

  local Vec fg2_fft[3];
  FOR_DIM
    fg2_fft[dim] = vec_pop (BHD->dc);

  /************************************************/
  /* rhob * FS*gab gb */
  /************************************************/

  /*
    This computes the k-space representation of the weighted force fft
    (f_ab  g_ab).  This  force  is the  precursor  to the  convolution
    kernel K_ab. Again the long-range part is treated separately.
  */
  {
    local Vec work = vec_pop (BHD->da);
    FOR_DIM
      {
        VecPointwiseMult (work, gab, fab_s[dim]);

        /* special treatment: Coulomb long */
        VecAXPY (work, -1.0, fab_l[dim]);

        MatMult (BHD->fft_mat, work, fg2_fft[dim]);
      }
    vec_push (BHD->da, &work);
  }

  {
    local Vec gb_fft = vec_pop (BHD->dc);

    /* The  convolution will  be  computed in  momentum  space, so  also
       compute fft (g_b): */
    MatMult (BHD->fft_mat, gb, gb_fft);

    /* This  effectively applies  K_ab  to g_b  in  k-space and  returns
       k-space du_a: */
    kapply (BHD, fg2_fft, gb_fft, cab_fft, dua_fft);

    vec_push (BHD->dc, &gb_fft);
  }

  FOR_DIM
    vec_push (BHD->dc, &fg2_fft[dim]);

  /* Transform the result of the convolution to real space du_a: */
  MatMultTranspose (BHD->fft_mat, dua_fft, dua);

  VecScale (dua, rhob * PD->beta/L/L/L);
}


/* w := x / max(y, thresh), essentially: */
static void safe_pointwise_divide (Vec w,        /* intent(out) */
                                   Vec x, Vec y) /* intent(in) */
{
  real pure f (real x, real y)
  {
    if (y < NORM_REG)
      return x / NORM_REG;
    else
      return x / y;
  }
  vec_map2 (w, f, x, y);
}


static double pure sinc (double x)
{
  if (unlikely (x == 0.0))
    return 1.0;
  else
    return sin (x) / x;
}


/*
  Output ω(k)  = sinc(kr), the  Fourier image of  ω(x) = δ(|x| -  r) /
  4πr².  As it appears that  time to compute sinc(k) is measurable, it
  may make sense  to precompute ω(k) for all  distinct r's. The origin
  is at the corner.
*/
static
void omega_intra (const ProblemData *PD, const DA dc, real r, Vec w_fft)
{
  const int *N = PD->N;         /* N[3] */
  const real L = PD->interval[1] - PD->interval[0];

  /* Get local portion of the k-grid */
  int x[3], n[3], i[3];
  DMDAGetCorners (dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***w_fft_;
  DMDAVecGetArray (dc, w_fft, &w_fft_);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);
          const real kr = (2.0 * M_PI * r / L) * sqrt (k2);

          /* Compute ω(k): */
          w_fft_[i[2]][i[1]][i[0]] = sinc (kr);
        }
  DMDAVecRestoreArray (dc, w_fft, &w_fft_);
}


/*
  Allocates and initializes  a matrix of intra-molecular correlations
  except of diagonal elements that are implicitly 1:
*/
void bgy3d_omega_fft_create (const State *BHD, int m, const Site solvent[m],
                             Vec omega_fft[m][m]) /* out, creates them */
{
  /* FIXME: m  x m  distance matrix does  not handle  equivalent sites
     well.  Diagonal zeros are never referenced: */
  real r[m][m];
  bgy3d_sites_dist_mat (m, solvent, r);

  /* Precompute omega[][]: */
  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)
      {
        omega_fft[j][i] = omega_fft[i][j] = bgy3d_vec_create (BHD->dc);
        omega_intra (BHD->PD, BHD->dc, r[i][j], omega_fft[i][j]);
        VecScale (omega_fft[j][i], BHD->PD->h[0] * BHD->PD->h[1] * BHD->PD->h[2]); /* historically */
      }

  /*
    One could  have allocated  them too and  fill with ones,  but this
    would waste  space. Insted a  NULL on the  diagonal implies ω  = 1
    identically:
  */
  for (int i = 0; i < m; i++)
    omega_fft[i][i] = NULL;
}


/*
  Output y(k)  is the  convolution of x(k)  with ω(x)  = δ(|x| -  r) /
  4πr².  This function takes momentum  space x(k) as input and returns
  momentum space y(k).  Historically,  this function operates of array
  of Vecs.   These arrays are  either of length  1 or 3,  in practice.
  All Vec's are complex here!

  NOTE: appears to work for in-place transformation with x == y.
*/
static void omega_apply (Vec w, int n, Vec x[n], Vec y[n])
{
  /* Set y(k) := w(k) * x(k). Note that w(k) is actually real: */
  complex pure f (complex w, complex x)
  {
    return w * x;
  }
  for (int i = 0; i < n; i++)
    vec_fft_map2 (y[i], f, w, x[i]);
}

/*
  Compute normalization function:

     c
    n  (x) = ∫ g  (y)  g  (x - y) dy
     ab         ac      bc

  where  one   of  the   convoluted  distributions  is   a  simplified
  "rigid-bond" distribution:  ω(x) = δ(|x|  - r) / 4πr².   This latter
  distribution is fully  characterized by the distance r  which is the
  only  input parameter. Another  distribution is  supplied as  a grid
  data.    NSSA   stays   for  "normalized   site-site   superposition
  approximation"  that  is   applied  to  triple  distributions.   See
  e.g. Eq. (4.106) on p. 75 of Jager thesis.

                                      c
  Thus, given g   and r   it returns n  .
               ac      bc             ab

  This function  takes momentum space g(k) as  input.  The convolution
  result is manipulated in  the real-space rep, unfortunately.  Though
  this  manipulation is plain  screening of  the too  small (negative)
  values.

  Vec nab is real, intent(out).

  Side effects: uses one complex temp Vec.

  FIXME: compare the code to Compute_dg_intra_ln().
 */
static void nssa_norm_intra (State *BHD, Vec gac_fft, Vec wbc_fft,
                             Vec nab)
{
  const real L = BHD->PD->interval[1] - BHD->PD->interval[0];

  local Vec nab_fft = vec_pop (BHD->dc);

  /* Set n(k) := ω(k) * g(k),  put result into Vec nab_fft.  Pass both
     gac_fft and nab_fft as arrays of length 1 to omega(): */
  omega_apply (wbc_fft, 1, &gac_fft, &nab_fft);

  /* Inverse FFT, n(k) -> n(x): */
  MatMultTranspose (BHD->fft_mat, nab_fft, nab);

  vec_push (BHD->dc, &nab_fft);

  VecScale (nab, 1./L/L/L);

  /* FIXME:  make  n(x) >  0,  because it  will  appear  later in  the
     denominator. Note that the same redundant precautions are made at
     other places too. */
  real pure f (real x)
  {
    const real thresh = 1.0e-8;
    if (x < thresh)
      return thresh;
    else
      return x;
  }
  /* Transfrom in-place: */
  vec_map1 (nab, f, nab); /* argument aliasing! */
}

/* Set the x[0, 0, 0] component of the Vec x to 0.0: */
static void set0 (Vec x)
{
  int ix[1] = {0};
  real y[1] = {0.0};

  /* Not collective. Should I call it from one worker only? */
  VecSetValues (x, 1, ix, y, INSERT_VALUES);
  VecAssemblyBegin (x);
  VecAssemblyEnd (x);
}

/*
  Compute  the  so-called  third,  "challenging",  term  of  the  NSSA
  expansion.  Given the input consisting of

         ~ b   c
    F  , g  , n  , and ω   represented by r
     ac   ac   ab       bc                 bc

  compute the divergence term to be added to

    u
     ab

  Long-range Coulomb is again  treated specially, thus the presence of
  the extra arguments: the long-range force fac_l and cab_fft.

  See line 3 of Eq. (4.118) p. 79 of Jager thesis.

  Side effects: uses one real and three complex work Vecs.
*/
static void Compute_dg_intra (State *BHD,
                              Vec fac[3], Vec fac_l[3],
                              Vec gac, Vec nab, Vec cac_fft,
                              Vec wbc_fft,
                              Vec dg_fft,
                              Vec dg, Vec dg_help)
{
  const ProblemData *PD = BHD->PD;
  const int *N = PD->N;         /* N[3] */

  const real h3 = PD->h[0] * PD->h[1] * PD->h[2];
  const real L = PD->interval[1] - PD->interval[0];
  const real scale = L / (2. * M_PI); /* siehe oben ... */

  local Vec fg2_fft[3];
  FOR_DIM
    fg2_fft[dim] = vec_pop (BHD->dc);

  /************************************************/
  /* Fa*ga ga*/
  /************************************************/

  /* fft(f1*gac) */
  {
    local Vec work = vec_pop (BHD->da);
    FOR_DIM
      {
        VecPointwiseMult (work, gac, fac[dim]);
        /* special treatment: Coulomb long */
        VecAXPY (work, -1.0, fac_l[dim]);

        MatMult (BHD->fft_mat, work, fg2_fft[dim]);
      }
    vec_push (BHD->da, &work);
  }

  /*
    Apply  omega() in-place,  note argument  aliasing.   Original code
    pretended that sinc(0) == 0. Is  it possible to claim that the k =
    0 component of weighted  forces vanishes?  Because then, formally,
    zero or one would not make a difference here.

    At least in HCl runs  the average of the short-range fac[] forces,
    which is the same as the  k = 0 component of its Fourier transform
    does  not vanish  for some  reason.  The  k =  0 component  of the
    long-range forces fac_l[] is zero.

    Now instead of redefining sinc(0) we  force the real part of k = 0
    component to 0:
  */
  FOR_DIM
    set0 (fg2_fft[dim]);
  omega_apply (wbc_fft, 3, fg2_fft, fg2_fft);

  /* int(..)/nab */
  /* Laplace^-1 * divergence */

  /*******************************************/
  /*    int/nab  */
  /*******************************************/

  /* Back  transformation of coulomb  part, divide  by nab  and forward
     transfromation */
  {
    local Vec work = vec_pop (BHD->da);
    FOR_DIM
      {
        MatMultTranspose (BHD->fft_mat, fg2_fft[dim], work);
        VecScale (work, 1./L/L/L);

        /* A   safer   version   of   VecPointwiseDivide   (BHD->v[dim],
           BHD->v[dim], nab);

           v[dim] := v[dim] / nab, essentially:
        */
        safe_pointwise_divide (work, work, nab);

        MatMult (BHD->fft_mat, work, fg2_fft[dim]);
      }
    vec_push (BHD->da, &work);
  }

  /* Get local portion of the k-grid */
  int x[3], n[3], i[3];
  DMDAGetCorners (BHD->dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  complex ***wbc_fft_, ***dg_fft_, ***cac_fft_;
  complex ***fg2_fft_[3];
  DMDAVecGetArray (BHD->dc, wbc_fft, &wbc_fft_);
  DMDAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  DMDAVecGetArray (BHD->dc, cac_fft, &cac_fft_);
  FOR_DIM
    DMDAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);

          /* Use fg2_fft: */
          complex div = 0.0;    /* complex */
          FOR_DIM
            div += ic[dim] * fg2_fft_[dim][i[2]][i[1]][i[0]];

          if (unlikely (k2 == 0))
            div *= 0.0;         /* 1/k2 is undefined */
          else
            div *= -I * ((h3 * scale) / k2);

          dg_fft_[i[2]][i[1]][i[0]] = div; /* complex */

          real sinc_kr;
          if (unlikely (k2 == 0))
            sinc_kr = 0.0; /* FIXME: sinc(0) == 1, but so the original */
          else
            sinc_kr = wbc_fft_[i[2]][i[1]][i[0]] / h3; /* lhs is real! */

          /* Overwrite fg2_fft: */
          FOR_DIM
            fg2_fft_[dim][i[2]][i[1]][i[0]] = ic[dim] / scale * (I * cac_fft_[i[2]][i[1]][i[0]]) * sinc_kr;
        }
  DMDAVecRestoreArray (BHD->dc, wbc_fft, &wbc_fft_);
  DMDAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  DMDAVecRestoreArray (BHD->dc, cac_fft, &cac_fft_);
  FOR_DIM
    DMDAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg_help);

  VecScale (dg_help, PD->beta/L/L/L);

  /* Back transformation  of coulomb part,  divide by nab  and forward
     transfromation */
  {
    local Vec work = vec_pop (BHD->da);
    FOR_DIM
      {
        MatMultTranspose (BHD->fft_mat, fg2_fft[dim], work);
        VecScale (work, 1./L/L/L);

        /* A safer version of VecPointwiseDivide (work, work, nab);

           work := work / nab, essentially:
        */
        safe_pointwise_divide (work, work, nab);

        MatMult (BHD->fft_mat, work, fg2_fft[dim]);
      }
    vec_push (BHD->da, &work);
  }

  DMDAVecGetArray (BHD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DMDAVecGetArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  /* loop over local portion of grid */
  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          /* FIXME: integer sum of squares will overflow for N >> 20000! */
          const int k2 = SQR (ic[2]) + SQR (ic[1]) + SQR (ic[0]);

          complex div = 0.0;    /* divergence */
          FOR_DIM
            div += ic[dim] * fg2_fft_[dim][i[2]][i[1]][i[0]];

          if (unlikely (k2 == 0))
            div *= 0.0;         /* 1/k2 is undefined */
          else
            div *= -I * ((h3 * scale) / k2);

          dg_fft_[i[2]][i[1]][i[0]] = div;
        }
  DMDAVecRestoreArray (BHD->dc, dg_fft, &dg_fft_);
  FOR_DIM
    DMDAVecRestoreArray (BHD->dc, fg2_fft[dim], &fg2_fft_[dim]);

  FOR_DIM
    vec_push (BHD->dc, &fg2_fft[dim]);

  MatMultTranspose (BHD->fft_mat, dg_fft, dg);

  VecScale (dg, PD->beta/L/L/L);

  VecAXPY (dg, 1.0, dg_help);
}


/*
  Compute   intramolecular  part.    The  first   thing  it   does  is
  transforming g(x) -> g(k). The real space representation g(x) is not
  otherwise used.

  Vec gac is intent(in).
  Vec dg is intent(out).

  Side   effects:  uses   one   complex  temp   Vec,   but  see   also
  nssa_norm_intra() which  uses one  more.  The two  can be  (and have
  been) aliased, but that is confusing.

  FIXME: compare the code to nssa_norm_intra().
*/
static void Compute_dg_intra_ln (State *BHD, Vec gac, Vec wbc_fft, Vec dg)
{
  local Vec gac_fft = vec_pop (BHD->dc);

  /* g(x) -> g(k): */
  MatMult (BHD->fft_mat, gac, gac_fft);

  /* Vec dg  is intent(out) here: */
  nssa_norm_intra (BHD, gac_fft, wbc_fft, dg);

  vec_push (BHD->dc, &gac_fft);

  /* -ln(g) */
  real pure f (real x)
  {
    // assert (x >= 1.0e-8);       /* See nssa_norm_intra() */
    return - log (x);
  }
  /* Transfrom in-place: */
  vec_map1 (dg, f, dg); /* argument aliasing! */

  /* Ensure normalization condition int(u)=0 */
  /* VecSum (dg, &k); */
  /* VecShift (dg, -k/N[0]/N[1]/N[2]); */
}


/*
  Approximate auxilary function

     c        ~ (2)              c
    Γ  (x) ~  g      = g  (x) / n  (x)
     ab        ab;c     ab       ab

  used to represent three-body distribution

                   c   b
    g    ~  g    Γ    Γ
     abc     bc   ab   ac

  in  normalized site-site superposition  approximation for  the mixed
  triplet distributions.   Of three sites the  two b, and  c belong to
  the  same rigid  species so  that their  pair distribution  is fully
  characterised by  the distance  r only.  See  e.g.  Eqs.   (4.102) -
  (4.106) on p. 75 of Jager thesis.

  Vec tab is real, intent(out).

  Side effects: by way of nssa_norm_intra() uses some work Vecs.
*/
static void nssa_gamma_cond (State *BHD, Vec gac_fft, Vec wbc_fft, Vec gab,
                             Vec tab) /* intent(out) */
{
  /* n(x) goes into Vec tab which is is intent(out) here: */
  nssa_norm_intra (BHD, gac_fft, wbc_fft, tab);

  /*
    t(x) =  g(x) / n(x)  (or rather t(x)  = g(x) / t(x)  with argument
    aliasing) avoiding small denominators.  Some of the commented code
    used VecPointwiseDivide() instead.
  */
  real pure f (real x, real y)
  {
    if (y < NORM_REG2)
      return x / NORM_REG2;
    else
      return x / y;
  }
  vec_map2 (tab, f, gab, tab); /* argument aliasing */
}

/*
  Compute  intramolecular  contribution due  to  b  /=  a.  So  called
  log-term, see e.g. Eqs.  (4.107), (4.109), and (4.114). Implies that
  the two  sites belong  to the  same species and  are separated  by a
  fixed distance behind intra-molecular correlation funciton wab_fft.

  To be called as in

  bgy3d_nssa_intra_log (BHD, g_fft[i], omega[i][j], g[j], du)

  The call sequence implemented by  this public wrapper also occurs in
  this  file at  a few  places.  It  appears that  one can  re-use the
  output of nssa_gamma_cond() though.
*/
void bgy3d_nssa_intra_log (State *BHD, Vec ga_fft, Vec wab_fft, Vec gb, Vec du)
{
  /*
    The  first  step  is  to  compute  normalization  functions  ĝ(x),
    Eqs.  (4.105),  (4.106)  and  (4.110).  This  already  involves  a
    convolution integral with the geometric factor δ(|x| - r) / 4πr².

    Vec du,  is intent (out) here and  is re-used as a  work vector to
    pass the result further.  Does one FFT^-1.
  */
  nssa_gamma_cond (BHD, ga_fft, wab_fft, gb, du);

  /*
    The next step is to  use the "conditional distribution" ĝ(x) again
    in a convolution  integral with the same geometric  factor δ(|x| -
    r) / 4πr², Eq. (4.114).

    Vec du  is both intent(in) in  first position and  intent (out) in
    the second.   Does one  FFT and one  FFT^-1.  Here it  is actually
    possible to use the same Vec as in- and output:
  */
  Compute_dg_intra_ln (BHD, du, wab_fft, du); /* aliasing! */
}


/* solve */
Vec BGY3d_solvent_solve (const ProblemData *PD, Vec g_ini)
{
  int m;                        /* number of solvent sites */
  const Site *solvent;          /* solvent[m] */

  assert(g_ini == PETSC_NULL);

  /* Get the number of solvent sites and their parameters: */
  bgy3d_solvent_get (&m, &solvent);

  bgy3d_solve_solvent (PD, m, solvent);

  return PETSC_NULL;

}

void bgy3d_solve_solvent (const ProblemData *PD, int m, const Site solvent[m])
{
  State *BHD = bgy3d_state_make (PD);

  /* Code used to be verbose: */
  bgy3d_problem_data_print (PD);

  int namecount=0;

  real du_norm[m][m];

  /* Original code used to print solvent params: */
  bgy3d_sites_show ("Solvent", m, solvent);

  PetscPrintf (PETSC_COMM_WORLD,
               "Solving BGY3d-M %d-site equation with Fourier ansatz...\n", m);

  /* allocation for local work vectors */
  Vec f[m][m][3], f_l[m][m][3]; /* full and long range pair force */

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      FOR_DIM
        {
          f[j][i][dim] = f[i][j][dim] = bgy3d_vec_create (BHD->da);
          f_l[j][i][dim] = f_l[i][j][dim] = bgy3d_vec_create (BHD->da);
        }

  /*
   The short-range solvent site-site  potentials for each pair (scaled
   by inverse  temperature beta) is  initially put into  the following
   array:
  */
  Vec u0[m][m];
  Vec c2[m][m];                 /* exp(- beta * LJ_repulsive(i, j) */

  bgy3d_vec_create2 (BHD->da, m, u0);
  bgy3d_vec_create2 (BHD->da, m, c2);

  /*
    Long-range Coulomb interaction for  solvent site pairs. So far the
    pairs differ only by a factor  q[i] * q[j]. Maybe we should rather
    store just one?  Here u2_fft are the Fourier transform of u2.  The
    same redundancy. Those are complex Vecs.
  */
  Vec u2[m][m], u2_fft[m][m];

  bgy3d_vec_create2 (BHD->da, m, u2);
  bgy3d_vec_create2 (BHD->dc, m, u2_fft); /* complex */

  /* Allocate more memory for fft */
  Vec gfg2_fft = bgy3d_vec_create (BHD->dc); /* complex */

  /* end of allocation */

  PetscPrintf (PETSC_COMM_WORLD, "Regularization of normalization: NORM_REG = %e\n", NORM_REG);
  PetscPrintf (PETSC_COMM_WORLD, "                                 NORM_REG2 = %e\n", NORM_REG2);

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Site specific  density.  Computed as a solvent  density rho times
     number of sites of that type in a solvent: */
  real rhos[m];
  for (int i = 0; i < m; i++)
    rhos[i] = PD->rho;

  /* Here dN is rho * dV, with dV being the weight of a grid point: */
  const real dN = PD->rho * PD->h[0] * PD->h[1] * PD->h[2];
  {
    const int N3 =  PD->N[0] * PD->N[1] * PD->N[2];
    PetscPrintf (PETSC_COMM_WORLD, "Number of solvent molecules is %f\n", N3 * dN);
  }

  /* Mixing parameter: */
  const real a0 = PD->lambda;

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  Vec g[m][m], du[m][m], t[m][m]; /* real, ij = ji */
  Vec g_fft[m][m];                /* complex, ij = ji */
  bgy3d_vec_create2 (BHD->da, m, g);
  bgy3d_vec_create2 (BHD->dc, m, g_fft); /* complex */
  bgy3d_vec_create2 (BHD->da, m, du);
  bgy3d_vec_create2 (BHD->da, m, t);

  Vec du_new = bgy3d_vec_create (BHD->da);
  Vec du_new2 = bgy3d_vec_create (BHD->da);
  Vec work = bgy3d_vec_create (BHD->da);

#ifdef L_BOUNDARY
  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem across iterations. See call to KSPSetInitialGuessNonzero()
    in bgy3d-poisson.c:
  */
  Vec x_lapl[m][m];             /* real */
  bgy3d_vec_create2 (BHD->da, m, x_lapl);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      VecSet (x_lapl[i][j], 0.0);
#endif

  local Vec omega[m][m];        /* diagonal will by NULL */
  bgy3d_omega_fft_create (BHD, m, solvent, omega); /* creates them */

  VecSet(du_new,0.0);

  for (real damp = damp_start; damp <= 1.0; damp += 0.1)
    {
      RecomputeInitialData (BHD, m, solvent,
                            f, f_l,
                            u0, c2,
                            u2, u2_fft,
                            (damp > 0.0 ? damp : 0.0), 1.0);
      PetscPrintf (PETSC_COMM_WORLD, "New lambda= %f\n", a0);

      /*
        Set initial guess, either here or by reading from file. At the
        end of the "dump" loop du[] is written to disk, so that in the
        next iteration we will read an updated version:
      */
      if (bgy3d_getopt_test ("--load-guess"))
        bgy3d_vec_read2 ("du%d%d.bin", m, du);
      else
        for (int i = 0; i < m; i++)
          for (int j = 0; j <= i; j++)
            VecSet (du[i][j], 0.0);

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          {
            bgy3d_impose_laplace_boundary (BHD, u0[i][j], x_lapl[i][j]);

            /* g = g0 * exp(-du) */
            compute_g (g[i][j], u0[i][j], du[i][j]);
          }

      /* Not sure if 0.0 as inital value is right. */
      real du_norm_old[m][m];
      for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
          du_norm_old[i][j] = 0.0;

      real a1 = a0;             /* loop-local variable */
  for (int iter = 0, mycount = 0, upwards = 0; iter < max_iter; iter++)
    {
      const int nth = 10;
      /*
        "a  = a1"  is taken  in  iteration 0,  10, 20,  etc.  "a1"  is
        modified during the loop.

        "a = a0" is taken in iterations 1-9, 11-19, etc.  "a0" remains
        unchanged during the loop.

        Note that in the first iteration a1 == a0.
      */

      /* Every nth  iteration, raise the mixing  coefficients just one
         time: */
      const real a = (iter % nth == 0) ? a1 : a0;

      /* Compute FFT of g[][] for all site pairs: */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          MatMult (BHD->fft_mat, g[i][j], g_fft[i][j]);

      /* For each site pair  ij (HH, HO, and OO) compute du  in g = g0
         exp(-du): */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          {
            /* Clear accumulator, terms will be added here: */
            VecSet (du_new, 0.0);

            /*
              Sum over  the third site k.   See the first  line of Eq.
              (4.118), p.  79 of  the Jager thesis.  This is so-called
              inter-molecular contribution. Also note that most of the
              pair quantities  are symmetric (I would say  all of them
              if I werent afraid of over-generalizations):

                Δu   = βρ  Σ  div (F   g  ) * g
                  ij     S  k       ik  ik     jk

              The site  density common for all sites  is replaced here
              by site-specific, but so far equal, density rho[k]. Note
              that in the 3-site code the sum over equivalent sites is
              treated  by scaling  the  density rho[k]  appropriately:
              e.g.  factor 2 for  H-site in H2O.  Equivalent sites are
              by  definition those  whose distributions  densities are
              equal.
            */
            for (int k = 0; k < m; k++)
              {
                Compute_dg_inter (BHD,
                                  f[i][k], f_l[i][k], g[i][k],
                                  g[j][k],
                                  u2_fft[i][k], rhos[k],
                                  gfg2_fft,
                                  du_new2);
                VecAXPY (du_new, 1.0, du_new2);
              }

            /* FIXME: 3-site code does not do this: */
            if (!bgy3d_getopt_test ("--no-hacks"))
              VecPointwiseMult (du_new, du_new, c2[i][j]);

            /*
              These are two sums  over k /= i and k /=  j. See lines 2
              and 3 of Eq. (4.118), p. 79 of Jager thesis:
            */
            for (int k = 0; k < m; k++)
              {
                /*
                  Line 2  of Eq.  (4.118),  p.  79 Jager  thesis. Here
                  sites k  and i belong  to the same  solvent species.
                  Intra-molecular  correlation is  fully  described by
                  the distance r[k][i] in the rigid solvent model.

                  Note that the site density rho[k] does not enter the
                  equation explicitly.   Instead the 3-site  code used
                  the literal  factor 2 when  incrementing accumulator
                  in a call to VecAXPY():
                */
                if (k != i)
                  {
                    /* Here t is intent(out): */
                    nssa_gamma_cond (BHD, g_fft[i][j], omega[k][i], g[j][k], t[j][k]);

                    /* Compute du_new2 term and add to the accumulator: */
                    Compute_dg_intra_ln (BHD, t[j][k], omega[k][i], du_new2);
                    VecAXPY (du_new, 1.0, du_new2);
                  }

                /*
                  Line  3 of Eq.   (4.118), p.   79 Jager  thesis. The
                  so-called "numerically challenging" term. Here sites
                  k  and  j  belong   to  the  same  solvent  species.
                  Intra-molecular  correlation is  fully  described by
                  the distance r[j][k] in the rigid solvent model:
                */
                if (k != j)
                  {               /* INTRA2 */
                    /*
                      We  need  two   distinct  temporaries  to  store
                      intermediates

                       ~ (2)        k
                       g      and  n  .
                        ik;j        ij

                      These two  qualify as k  != j in the  context of
                      this summation:
                    */
                    assert (t[i][k] != t[i][j]);

                    /*
                      NOTE:  without  the  next conditional  the  next
                      computation would be redundant computation for j
                      == i.  See the call to nssa_gamma_cond() above:
                    */
                    if (i != j)
                      nssa_gamma_cond (BHD, g_fft[i][j], omega[j][k], g[i][k], t[i][k]);
                    nssa_norm_intra (BHD, g_fft[i][k], omega[j][k], t[i][j]);
                    Compute_dg_intra (BHD,
                                      f[i][k], f_l[i][k], t[i][k],
                                      t[i][j],
                                      u2_fft[i][k], omega[j][k],
                                      gfg2_fft,
                                      du_new2, work);
                    VecAXPY (du_new, 1.0, du_new2);
                  }
              }

            /* Long-range Coulomb, scaled by inverse temperature: */
            VecAXPY (du_new, PD->beta, u2[i][j]);

            /*
              Add  an  effective  Coulomb  field of  boundary  surface
              charge  to  make  the  "potential" du  at  the  boundary
              vanish:
            */
            if (iter >= 0)
              bgy3d_impose_laplace_boundary (BHD, du_new, x_lapl[i][j]);

            /*
              Mix du and du_new with a fixed ratio "a":

                du' = a * du_new + (1 - a) * du
                norm = |du_new - du|

              last arg is a temp
            */
            du_norm[i][j] = bgy3d_vec_mix (du[i][j], du_new, a, work);
            du_norm[j][i] = du_norm[i][j];
          }

      /*
        Now  that du[]  has been  computed using  g[] of  the previous
        iteration one can safely update  g[].  Compute g := exp[-(u0 +
        du)], with a sanity check:
      */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          compute_g (g[i][j], u0[i][j], du[i][j]);

      /* Fancy step  size control. FIXME:  weired logic. Code  used to
         check if *any* of the norms went up: */
      bool up;
      {
        real diff[m][m];
        for (int i = 0; i < m; i++)
          for (int j = 0; j < m; j++)
            diff[i][j] = du_norm[i][j] - du_norm_old[i][j];
        up = maxval (m * m, (real*) diff) > 0.0;
      }

      /* That was the only place comparing to du_norm_old[]: */
      for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) /* Yes, full range! */
          du_norm_old[i][j] = du_norm[i][j];

      mycount++;

      if (iter % nth != 1 && up) /* not in the nth + 1 iteration */
        upwards = 1;
      else if (iter > 2 * nth && iter % nth == 1 && upwards == 0 && up)
        {
          /* In the  nth + 1 iteration,  if the norm  went up decrease
             the mixing: */
          a1 = MAX (a1 / 2.0, a0);
          mycount = 0;
        }
      else
        upwards = 0;

      /* Scale the coefficient  "a1" up by a factor,  but make sure it
         is not above 1.0. Reset mycount. */
      if (mycount > 2 * nth)
        {
          a1 = MIN (a1 * 2.0, 1.0);
          mycount = 0;
        }
      /* otherwise leave "a1" and "mycount" unchanged */

      PetscPrintf (PETSC_COMM_WORLD, "%03d ", iter + 1);
      PetscPrintf (PETSC_COMM_WORLD, "a=%f ", a);

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          PetscPrintf (PETSC_COMM_WORLD, "%s-%s=%e ",
                       solvent[i].name, solvent[j].name, du_norm[i][j]);

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          PetscPrintf (PETSC_COMM_WORLD,
                       "h(%s-%s)=% f ",
                       solvent[i].name,
                       solvent[j].name,
                       dN * vec_hole (g[i][j]));

      PetscPrintf (PETSC_COMM_WORLD, "count=%3d upwards=%1d", mycount, upwards);
      PetscPrintf (PETSC_COMM_WORLD, "\n");

      /* Exit  when  any  of  du[]   does  not  change  by  more  than
         norm_tol: */
      if (maxval (m * m, (real*) du_norm) <= norm_tol)
        {
          PetscPrintf (PETSC_COMM_WORLD,
                       "norm %e <= %e (norm-tol) in iteration %d < %d (max-iter)\n",
                       maxval (m * m, (real*) du_norm), norm_tol, iter + 1, max_iter);
          break;
        }
    } /* for (iter = ... ) */

  /* FIXME: Debug  output from every iteration  with different overall
     scale factors damp.  Remove when no more needed. */
  {
    char fmt[20];
    snprintf (fmt, sizeof fmt, "vec%%d%%d-%d.m", namecount++);
    bgy3d_vec_save_ascii2 (fmt, m, g);
  }

  /* Save du[][] to binary files: */
  if (bgy3d_getopt_test ("--save-guess"))
    bgy3d_vec_save2 ("du%d%d.bin", m, du);

  /* Save g2[][] to binary files: */
  bgy3d_vec_save2 ("g%d%d.bin", m, g);

    } /* for (dump = ... ) */

  /* deallocation of work vectors */
  bgy3d_vec_destroy2 (m, g);
  bgy3d_vec_destroy2 (m, g_fft);
  bgy3d_vec_destroy2 (m, t);
  bgy3d_vec_destroy2 (m, du);
  bgy3d_vec_destroy2 (m, x_lapl);

  bgy3d_vec_destroy (&du_new);
  bgy3d_vec_destroy (&du_new2);
  bgy3d_vec_destroy (&work);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < i; j++)
      bgy3d_vec_destroy (&omega[i][j]);

  bgy3d_vec_destroy2 (m, u2);
  bgy3d_vec_destroy2 (m, u2_fft);
  bgy3d_vec_destroy2 (m, u0);
  bgy3d_vec_destroy2 (m, c2);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        bgy3d_vec_destroy1 (3, f[i][j]);
        bgy3d_vec_destroy1 (3, f_l[i][j]);
      }

  bgy3d_vec_destroy (&gfg2_fft);

  bgy3d_state_destroy (BHD);
}


/* solve with product ansatz g=g0*dg */
Vec BGY3d_solvent_solve_h2o (const ProblemData *PD, Vec g_ini)
{
  State *BHD = bgy3d_state_make (PD);

  Vec du_new, du_new2, work;
  int namecount=0;

  int m;                        /* number of solvent sites */
  const Site *solvent;          /* solvent[m] */

  /* Get the number of solvent sites and their parameters: */
  bgy3d_solvent_get (&m, &solvent);

  /* Original code used to print solvent params: */
  bgy3d_sites_show ("Solvent", m, solvent);

  /* Original code only used arrays like F[2][2] */
  assert (m == 2);

  real du_norm[m][m];

  assert(g_ini == PETSC_NULL);

  PetscPrintf(PETSC_COMM_WORLD, "Solving BGY3dM (H2O) equation with Fourier ansatz...\n");

  /* allocation for local work vectors */
  Vec f[m][m][3], f_l[m][m][3];
  Vec u0[m][m], c2[m][m];
  Vec u2[m][m], u2_fft[m][m];

  bgy3d_vec_create2 (BHD->da, m, u2);
  bgy3d_vec_create2 (BHD->dc, m, u2_fft); /* complex */
  bgy3d_vec_create2 (BHD->da, m, c2);
  bgy3d_vec_create2 (BHD->da, m, u0);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      FOR_DIM
        {
          f[j][i][dim] = f[i][j][dim] = bgy3d_vec_create (BHD->da);
          f_l[j][i][dim] = f_l[i][j][dim] = bgy3d_vec_create (BHD->da);
        }

  /* Allocate more memory for fft */
  Vec gfg2_fft = bgy3d_vec_create (BHD->dc); /* complex */

  /* end of allocation */

  if (r_HH < 0)
    {
      PetscPrintf (PETSC_COMM_WORLD, "Solvent not a 3-Site model!\n");
      exit(1);
    }

  /*
   * Extract BGY3d specific things from supplied input:
   */

  /* Site specific  density.  Computed as a solvent  density rho times
     number of sites of that type in a solvent: */
  real rhos[m];
  for (int i = 0; i < m; i++)
    rhos[i] = PD->rho;
  rhos[0] = 2.0 * rhos[0];

  /* Mixing parameter: */
  const real a0 = PD->lambda;

  /* Initial damping factor: */
  const real damp_start = PD->damp;

  /* Number of total iterations */
  const int max_iter = PD->max_iter;

  /* norm_tol for convergence test */
  const real norm_tol = PD->norm_tol;

  Vec g[m][m], du[m][m], t[m][m]; /* real, ij = ji */
  Vec g_fft[m][m];                /* complex, ij = ji */
  bgy3d_vec_create2 (BHD->da, m, g);
  bgy3d_vec_create2 (BHD->dc, m, g_fft); /* complex */
  bgy3d_vec_create2 (BHD->da, m, du);
  bgy3d_vec_create2 (BHD->da, m, t);

  du_new = bgy3d_vec_create (BHD->da);
  du_new2 = bgy3d_vec_create (BHD->da);
  work = bgy3d_vec_create (BHD->da);

#ifdef L_BOUNDARY
  /*
    These  will be  used to  store solutions  of the  Laplace boundary
    problem across iterations. See call to KSPSetInitialGuessNonzero()
    in bgy3d-poisson.c:
  */
  Vec x_lapl[2][2];             /* real */
  bgy3d_vec_create2 (BHD->da, 2, x_lapl);

  for (int i = 0; i < 2; i++)
    for (int j = 0; j <= i; j++)
      VecSet (x_lapl[i][j], 0.0);
#endif

  Vec omega[2][2];
  omega[0][1] = bgy3d_vec_create (BHD->dc);
  omega[0][0] = bgy3d_vec_create (BHD->dc);
  omega[1][0] = omega[0][1];
  omega[1][1] = NULL;

  omega_intra (BHD->PD, BHD->dc, r_HO, omega[0][1]);
  omega_intra (BHD->PD, BHD->dc, r_HH, omega[0][0]);

  VecScale (omega[0][1], BHD->PD->h[0] * BHD->PD->h[1] * BHD->PD->h[2]); /* historically */
  VecScale (omega[0][0], BHD->PD->h[0] * BHD->PD->h[1] * BHD->PD->h[2]); /* historically */

  /* Set initial guess, either here or by reading from file: */
  if (bgy3d_getopt_test ("--load-guess"))
    bgy3d_vec_read2 ("du%d%d.bin", m, du);
  else
    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i; j++)
        VecSet (du[i][j], 0.0);


  VecSet(du_new,0.0);

  for (real damp = damp_start; damp <= damp_start; damp += 0.1)
    {
      RecomputeInitialData (BHD, m, solvent,
                            f, f_l,
                            u0, c2,
                            u2, u2_fft,
                            (damp > 0.0 ? damp : 0.0), 1.0);
      PetscPrintf (PETSC_COMM_WORLD, "New lambda= %f\n", a0);

      bgy3d_impose_laplace_boundary (BHD, u0[0][0], x_lapl[0][0]);
      bgy3d_impose_laplace_boundary (BHD, u0[1][1], x_lapl[1][1]);
      bgy3d_impose_laplace_boundary (BHD, u0[0][1], x_lapl[0][1]);

      /* g=g0*exp(-du) */
      compute_g (g[0][1], u0[0][1], du[0][1]);
      compute_g (g[0][0], u0[0][0], du[0][0]);
      compute_g (g[1][1], u0[1][1], du[1][1]);

      /* Not sure if 0.0 as inital value is right. */
      real du_norm_old[m][m];
      for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
          du_norm_old[i][j] = 0.0;

      real a1 = a0;             /* loop-local variable */
  for (int iter = 0, mycount = 0, upwards = 0; iter < max_iter; iter++)
    {
      /* Every nth iteration, starting with iter == 0: */
      const bool nth = !(iter % 10);

      /*
        "a  = a1"  is taken  in  iteration 0,  10, 20,  etc.  "a1"  is
        modified during the loop.

        "a = a0" is taken in iterations 1-9, 11-19, etc.  "a0" remains
        unchanged during the loop.

        Note that in the first iteration a1 == a0.
      */
      const real a = nth? a1 : a0;

      /* Compute FFT of g[][] for all site pairs: */
      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          MatMult (BHD->fft_mat, g[i][j], g_fft[i][j]);

      /* f=integral(g) */
      if (1)                    /* kflg was set when with -pair */
        {
          /* g_OH */
          VecSet (du_new, 0.0);
          Compute_dg_inter (BHD,
                            f[1][1], f_l[1][1], g[1][1],
                            g[0][1],
                            u2_fft[1][1], rhos[1],
                            gfg2_fft, du_new2);
          VecAXPY (du_new, 1.0, du_new2);

          Compute_dg_inter (BHD,
                            f[0][1], f_l[0][1], g[0][1],
                            g[0][0],
                            u2_fft[0][1], rhos[0],
                            gfg2_fft, du_new2);
          VecAXPY (du_new, 1.0, du_new2);

          nssa_gamma_cond (BHD, g_fft[0][1], omega[0][1], g[0][0], t[0][0]);
          Compute_dg_intra_ln (BHD, t[0][0], omega[0][1], du_new2); /* t is intent(in) */
          VecAXPY(du_new, 2.0, du_new2);

          /* tO = gHO/int(gHO wHH) */
          nssa_gamma_cond (BHD, g_fft[0][1], omega[0][0], g[0][1], t[1][1]);
          nssa_norm_intra (BHD, g_fft[0][1], omega[0][0], t[0][1]);
          Compute_dg_intra (BHD, f[0][1], f_l[0][1],
                            t[1][1], t[0][1],
                            u2_fft[0][1], omega[0][0], gfg2_fft, du_new2, work);
          VecAXPY(du_new, 1.0, du_new2);

          nssa_gamma_cond (BHD, g_fft[0][1], omega[0][1], g[1][1], t[1][1]);
          nssa_norm_intra (BHD, g_fft[1][1], omega[0][1], t[0][1]);
          Compute_dg_intra (BHD, f[1][1], f_l[1][1],
                            t[1][1], t[0][1],
                            u2_fft[1][1], omega[0][1], gfg2_fft, du_new2, work);
          VecAXPY(du_new, 1.0, du_new2);

          VecAXPY(du_new, PD->beta, u2[0][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, du_new, x_lapl[0][1]);

          du_norm[0][1] = bgy3d_vec_mix (du[0][1], du_new, a, work);

          /* g_H */
          VecSet (du_new, 0.0);
          Compute_dg_inter (BHD,
                            f[0][1], f_l[0][1], g[0][1],
                            g[0][1],
                            u2_fft[0][1], rhos[1],
                            gfg2_fft, du_new2);
          VecAXPY (du_new, 1.0, du_new2);

          Compute_dg_inter (BHD,
                            f[0][0], f_l[0][0], g[0][0],
                            g[0][0],
                            u2_fft[0][0], rhos[0],
                            gfg2_fft, du_new2);
          VecAXPY (du_new, 1.0, du_new2);

          nssa_gamma_cond (BHD, g_fft[0][0], omega[0][0], g[0][0], t[0][0]);
          Compute_dg_intra_ln (BHD, t[0][0], omega[0][0], du_new2); /* t is intent(in) */
          VecAXPY(du_new, 1.0, du_new2);

          nssa_gamma_cond (BHD, g_fft[0][0], omega[0][1], g[0][1], t[0][1]);
          Compute_dg_intra_ln (BHD, t[0][1], omega[0][1], du_new2); /* t is intent(in) */
          VecAXPY(du_new, 1.0, du_new2);

          /* tO = gH/int(gH wHH) */
          nssa_gamma_cond (BHD, g_fft[0][0], omega[0][0], g[0][0], t[1][1]);
          nssa_norm_intra (BHD, g_fft[0][0], omega[0][0], t[0][0]);
          Compute_dg_intra (BHD, f[0][0], f_l[0][0],
                            t[1][1], t[0][0],
                            u2_fft[0][0], omega[0][0], gfg2_fft, du_new2, work);
          VecAXPY(du_new, 1.0, du_new2);

          nssa_gamma_cond (BHD, g_fft[0][0], omega[0][1], g[0][1], t[0][1]);
          nssa_norm_intra (BHD, g_fft[0][1], omega[0][1], t[0][0]);
          Compute_dg_intra (BHD, f[0][1], f_l[0][1],
                            t[0][1], t[0][0],
                            u2_fft[0][1], omega[0][1], gfg2_fft, du_new2, work);
          VecAXPY(du_new, 1.0, du_new2);

          VecAXPY(du_new, PD->beta, u2[0][0]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, du_new, x_lapl[0][0]);

          du_norm[0][0] = bgy3d_vec_mix (du[0][0], du_new, a, work);

          /* g_O */
          VecSet (du_new, 0.0);
          Compute_dg_inter (BHD,
                            f[0][1], f_l[0][1], g[0][1],
                            g[0][1],
                            u2_fft[0][1], rhos[0],
                            gfg2_fft, du_new2);
          VecAXPY (du_new, 1.0, du_new2);

          Compute_dg_inter (BHD,
                            f[1][1], f_l[1][1], g[1][1],
                            g[1][1],
                            u2_fft[1][1], rhos[1],
                            gfg2_fft, du_new2);
          VecAXPY (du_new, 1.0, du_new2);


          nssa_gamma_cond (BHD, g_fft[1][1], omega[0][1], g[0][1], t[0][1]);
          Compute_dg_intra_ln (BHD, t[0][1], omega[0][1], du_new2); /* t is intent(in) */
          VecAXPY(du_new, 2.0, du_new2);

          nssa_gamma_cond (BHD, g_fft[1][1], omega[0][1], g[0][1], t[0][1]);
          nssa_norm_intra (BHD, g_fft[0][1], omega[0][1], t[1][1]);
          Compute_dg_intra (BHD, f[0][1], f_l[0][1],
                            t[0][1], t[1][1],
                            u2_fft[0][1], omega[0][1], gfg2_fft, du_new2, work);
          VecAXPY(du_new, 2.0, du_new2);

          VecAXPY(du_new, PD->beta, u2[1][1]);

          if (iter >= 0)
            bgy3d_impose_laplace_boundary (BHD, du_new, x_lapl[1][1]);

          du_norm[1][1] = bgy3d_vec_mix (du[1][1], du_new, a, work);

          /* ende: */
          compute_g (g[0][1], u0[0][1], du[0][1]);
          compute_g (g[0][0], u0[0][0], du[0][0]);
          compute_g (g[1][1], u0[1][1], du[1][1]);
        } /* of if (1) */

      /* (fancy) step size control */
      mycount++;
      if (((iter-1)%10) &&
          (du_norm_old[0][0] < du_norm[0][0] ||
           du_norm_old[1][1] < du_norm[1][1] ||
           du_norm_old[0][1] < du_norm[0][1]))
        {
          upwards = 1;
        }
      else if (iter>20 && !((iter-1)%10) && upwards==0 &&
               (du_norm_old[0][0] < du_norm[0][0] ||
                du_norm_old[1][1] < du_norm[1][1] ||
                du_norm_old[0][1] < du_norm[0][1]))
        {
          a1 /= 2.0;
          if (a1 < a0)
            a1 = a0;
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

      for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++) /* Yes, full range! */
          du_norm_old[i][j] = du_norm[i][j];

      PetscPrintf (PETSC_COMM_WORLD, "%03d ", iter + 1);
      PetscPrintf (PETSC_COMM_WORLD, "a=%f ", a);

      for (int i = 0; i < m; i++)
        for (int j = 0; j <= i; j++)
          PetscPrintf (PETSC_COMM_WORLD, "%s-%s=%e ",
                       solvent[i].name, solvent[j].name, du_norm[i][j]);

      PetscPrintf (PETSC_COMM_WORLD, "count=%3d upwards=%1d", mycount, upwards);
      PetscPrintf (PETSC_COMM_WORLD, "\n");

      if (du_norm[0][0] <= norm_tol &&
          du_norm[1][1] <= norm_tol &&
          du_norm[0][1] <= norm_tol)
        break;
    } /* for (iter = ...) */

  /* FIXME: Debug  output from every iteration  with different overall
     scale factors damp.  Remove when no more needed. */
  {
    char fmt[20];
    snprintf (fmt, sizeof fmt, "vec%%d%%d-%d.m", namecount++);
    bgy3d_vec_save_ascii2 (fmt, m, g);
  }

  /* Save du[][] to binary files: */
  if (bgy3d_getopt_test ("--save-guess"))
    bgy3d_vec_save2 ("du%d%d.bin", 2, du);

  /* Save g2[][] to binary files: */
  bgy3d_vec_save2 ("g%d%d.bin", 2, g);

    } /* for (damp = ...) */

  bgy3d_vec_destroy2 (m, g);
  bgy3d_vec_destroy2 (m, g_fft);
  bgy3d_vec_destroy2 (m, t);
  bgy3d_vec_destroy2 (m, du);
  bgy3d_vec_destroy2 (m, x_lapl);

  bgy3d_vec_destroy (&du_new);
  bgy3d_vec_destroy (&du_new2);
  bgy3d_vec_destroy (&work);

  bgy3d_vec_destroy (&omega[0][1]);
  bgy3d_vec_destroy (&omega[0][0]);

  bgy3d_vec_destroy2 (m, u2);
  bgy3d_vec_destroy2 (m, u2_fft);
  bgy3d_vec_destroy2 (m, u0);
  bgy3d_vec_destroy2 (m, c2);

  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        bgy3d_vec_destroy1 (3, f[i][j]);
        bgy3d_vec_destroy1 (3, f_l[i][j]);
      }

  bgy3d_vec_destroy (&gfg2_fft);

  bgy3d_state_destroy (BHD);


  return PETSC_NULL;
}

