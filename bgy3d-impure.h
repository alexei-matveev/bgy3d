/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013, 2014 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/*
  A Restart* token to be  passed there an back for resuming iterations
  without  disk access.  There is  no explicit  constructor  for these
  objects     but    they     are    (optionally)     returned    from
  bgy3d_solute_solve().
*/
typedef struct Restart Restart;

/*
  This function is the main entry  point for the BGY3dM equation for a
  m-site solvent and an arbitrary solute.  The vectors in

  Vec g[m], intent(out)

  are  initialzed as  global distributed  arrays and  filled  with the
  solvent site  distributions. It is the responsibility  of the caller
  to destroy them when no more needed.
*/
void bgy3d_solute_solve (const ProblemData *PD,
                         int m, const Site solvent[m],
                         int n, const Site solute[n],
                         void (*density)(int k, const real x[k][3], real rho[k]),
                         Vec g[m],          /* out */
                         Context **medium,  /* out, optional */
                         Restart **restart); /* inout, optional */

void bgy3d_restart_destroy (Restart *restart);

Vec BGY3d_solute_solve (const ProblemData *PD, Vec g_ini);
