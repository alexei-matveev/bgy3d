/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-solutes.h"      /* Site */
#include "bgy3d-solvents.h"     /* FIXME: many #defines */

/*
 These  are the  two solvent  sites.   Coordinates will  not be  used.
 Respective parameters are #defined  in bgy3d-solvents.h.  Also do not
 take the names of the sites literally.  The same structure is used to
 represent all (2-site) solvents, such as HCl.

 FIXME: at many places it is  assumed that the number of solvent sites
 is exactly two:
*/
static const Site solvent[] =
  {{"h", {0.0, 0.0, 0.0}, sH, eH, qH}, /* dont use sH, eH, qH directly */
   {"o", {0.0, 0.0, 0.0}, sO, eO, qO}}; /* same for sO, eO, qO */

#undef sH
#undef eH
#undef qH
#undef sO
#undef eO
#undef qO


/*
  Get solvent sites.  Functions that do the real work operate on array
  of sites.  As  of now the returned sites[] do not  need to be freed.
*/
void bgy3d_solvent_get (/* const char *name */ int *n, const Site **sites)
{
  /* Return solvent sites and their count: */
  *n = sizeof (solvent) / sizeof (Site);
  *sites = solvent;
}

/*
  Prints a table like this:

Solvent:
#       site    ε               σ                q
1       h       0.039710        2.735000         0.200000
2       o       0.514340        3.353000        -0.200000
*/
void bgy3d_sites_show (const char *name, int m, const Site sites[m])
{
  PetscPrintf (PETSC_COMM_WORLD, "%s:\n", name);
  PetscPrintf (PETSC_COMM_WORLD, "#\tsite\tε       \tσ       \t q\n"); /* unicode here! */
  for (int i = 0; i < m; i++)
    PetscPrintf (PETSC_COMM_WORLD,
                 "%d\t%-4s\t%f\t%f\t% f\n",
                 i + 1,
                 sites[i].name,
                 sites[i].epsilon,
                 sites[i].sigma,
                 sites[i].charge);
}

