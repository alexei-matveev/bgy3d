/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

/*
  A function  that takes a  context, an input  Vec x and  computes the
  residual Vec r.  An example is the "error  vector" of the non-linear
  equation is  the difference  between the input  and the output  of a
  "fixpoint" iteration as a function of input:

    r = x    - x
         out    in
*/
typedef void (*Function) (void *ctx, Vec x, Vec r);

/* Solver for  non-linear equations, either Newton  or fixpoint Picard
   iterations: */
typedef void (*Solver) (const ProblemData *PD, void *ctx, Function F, Vec x);

/*
  A few solvers for F(x) = 0 taking a Function, its execution context,
  initial guess, and some  user input eventually affecting convergence
  criteria:
*/
void bgy3d_snes_default (const ProblemData *PD, void *ctx, Function F, Vec x);
void bgy3d_snes_newton (const ProblemData *PD, void *ctx, Function F, Vec x);
void bgy3d_snes_picard (const ProblemData *PD, void *ctx, Function F, Vec x);
void bgy3d_snes_jager (const ProblemData *PD, void *ctx, Function F, Vec x);
