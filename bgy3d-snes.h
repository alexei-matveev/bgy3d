/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/*
  A function  that takes a  context, an input  Vec x and  computes the
  residual Vec r.  An example is the "error  vector" of the non-linear
  equation is  the difference  between the input  and the output  of a
  "fixpoint" iteration as a function of input:

    r = x    - x
         out    in
*/
typedef void (*ArrFunc1) (void *ctx, int n, const real x[n],
                            real r[n]); /* out */

typedef void (*ArrFunc2) (void *ctx, int n, const real x[n],  const real y[n],
                            real r[n]); /* out */

typedef void (*VecFunc1) (void *ctx, Vec x, Vec r /* out */);

/* A function to apply Jacobian: r = J(x) * dx */
typedef void (*VecFunc2) (void *ctx, Vec x, Vec dx, Vec r /* out */);

/*
  Solvers for  non-linear equations, either Newton  or fixpoint Picard
  iterations.

  A solver for  an untyped (array) form ArrFunc1.   Finds an x_[] such
  that  dx_[]  as  returned  by  ArrFunc1  f  (ctx,  n,  x_,  dx_)  is
  zero. ArrFunc2 df, when not NULL applies the Jacobian.

  FIXME: also adapt interface declaration in ./snes.f90!
*/
void rism_snes (void *ctx, ArrFunc1 f, ArrFunc2 df, int n, real x_[n]);

/* For solving linear equation F(x) = b iteratively. */
void bgy3d_krylov (void *ctx, VecFunc1 F, Vec b, Vec x);

/* Solves for f(x) = b  iteratively. Has to be consistent with Fortran
   declarations in snes.f90: */
void rism_krylov (void *ctx, ArrFunc1 f, int n, real b_[n], real x_[n]);

/*
  A  few solvers  for  F(x) =  0  taking a  VecFunc1, its  execution
  context,  initial guess,  and some  user input  eventually affecting
  convergence criteria:
*/
void bgy3d_snes_default (const ProblemData *PD, void *ctx,
                         VecFunc1 F, VecFunc2 dF, Vec x);

void bgy3d_snes_newton (const ProblemData *PD, void *ctx,
                        VecFunc1 F, VecFunc2 dF, Vec x);

void bgy3d_snes_picard (const ProblemData *PD, void *ctx,
                        VecFunc1 F, VecFunc2 dF, Vec x);

void bgy3d_snes_jager (const ProblemData *PD, void *ctx,
                       VecFunc1 F, VecFunc2 dF, Vec x);

void bgy3d_snes_trial (const ProblemData *PD, void *ctx,
                       VecFunc1 F, VecFunc2 dF, Vec x);
