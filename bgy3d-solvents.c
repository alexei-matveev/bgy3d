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
