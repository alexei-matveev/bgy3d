/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/*
  A function  that takes a  context, an input  Vec x and  computes the
  residual Vec r.  An example is the "error  vector" of the non-linear
  equation is  the difference  between the input  and the output  of a
  "fixpoint" iteration as a function of input:

    r = x    - x
         out    in
*/
typedef void (*ArrayFunc) (void *ctx, int n,  const real x[n],  real r[n]);
typedef void (*VectorFunc) (void *ctx, /* const */ Vec x, /* out */ Vec r);

/*
  Solvers for  non-linear equations, either Newton  or fixpoint Picard
  iterations:

  typedef void (*VectorSolver) (const ProblemData *PD, void *ctx, VectorFunc F, Vec x);
  typedef void (*ArraySolver) (void *ctx, ArrayFunc f, int n, real x[n]);
*/

/*
  A solver for  an untyped (array) form ArrayFunc.  Finds an x_[] such
  that dx_[] as returned by ArrayFunc f (ctx, n, x_, dx_) is zero.
*/
void rism_snes (void *ctx, ArrayFunc f, int n, real x_[n]);

/*
  A  few solvers  for  F(x) =  0  taking a  VectorFunc, its  execution
  context,  initial guess,  and some  user input  eventually affecting
  convergence criteria:
*/
void bgy3d_snes_default (const ProblemData *PD, void *ctx, VectorFunc F, Vec x);
void bgy3d_snes_newton (const ProblemData *PD, void *ctx, VectorFunc F, Vec x);
void bgy3d_snes_picard (const ProblemData *PD, void *ctx, VectorFunc F, Vec x);
void bgy3d_snes_jager (const ProblemData *PD, void *ctx, VectorFunc F, Vec x);
void bgy3d_snes_trial (const ProblemData *PD, void *ctx, VectorFunc F, Vec x);
