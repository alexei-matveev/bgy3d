/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*
  Copyright (c) 2007 Lukas Jager
  Copyright (c) 2013 Alexei Matveev
  Copyright (c) 2013 Bo Li
*/

/*==========================================*/
/* Hinweise zur Faltung: */
/* g(x) muss immer um 0 herum bekannt sein. Also sollten auch */
/* die x_M um 0 herum liegen. */
/*===============================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3d-solutes.h"      /* FIXME: struct Site in prototypes */
#include "bgy3d-pure.h"         /* BGY3d_solvent_solve(), ... */
#include "bgy3d-potential.h"    /* Context */
#include "bgy3d-impure.h"       /* BGY3d_solute_solve() */
#include "bgy3d-vec.h"          /* bgy3d_vec_*() */
#include "hnc3d.h"              /* hnc3d_solute_solve(), ... */

#ifdef WITH_GUILE
#include "bgy3d-guile.h"        /* bgy3d_guile_main */
#endif

#ifdef WITH_GUILE
int main (int argc, char **argv)
{
  bgy3d_guile_main (argc, argv);
}
#else
static char helptext[] = "Solving BGY3d equation.\n";

typedef Vec Solver (const ProblemData *PD, Vec g_ini);

int main (int argc, char **argv)
{
  real mpi_start, mpi_stop;
  Solver *solver = NULL;

  verbosity = 0;                /* global var */

  /* Petsc  or MPI  may  choose to  rewrite  the command  line, do  it
     early: */
  PetscInitialize (&argc, &argv, (char*)0, helptext);

  /* Make Petsc abort when it encounters an error: */
  PetscPushErrorHandler (PetscAbortErrorHandler, NULL);

  /*
    Read the  command line options.  Petsc insists on keys  having the
    leading dash, so  keep them for the moment.   Set global verbosity
    early  enough. This is  the only  short option!   Use long-options
    prefixed by "--" as the usual convention.
  */
  bgy3d_getopt_int ("-v", &verbosity);

  /* This calls bgy3d_getopt_*() a few more times: */
  ProblemData PD = bgy3d_problem_data ();

  /* computation time measurement start point */
  MPI_Barrier (PETSC_COMM_WORLD);
  mpi_start = MPI_Wtime();

  PetscViewerSetFormat(PETSC_VIEWER_STDERR_WORLD, PETSC_VIEWER_ASCII_MATLAB);
  PetscViewerSetFormat(PETSC_VIEWER_STDERR_SELF, PETSC_VIEWER_ASCII_MATLAB);
  /*==================================*/

  {
    int np;
    MPI_Comm_size (PETSC_COMM_WORLD, &np);
    PetscPrintf (PETSC_COMM_WORLD, "NP= %d\n", np);
  }

  /* FIXME:  some solvers do  this, some  dont. Sometimes  this output
     appears twice: */
  bgy3d_problem_data_print (&PD);

  //PetscPrintf(PETSC_COMM_WORLD, "\tATTENTION: Factor 2 is included!!! But why???\n");

  /* if(PD.id==1) */
  /*     start_debugger(); */
  /* sleep(5); */

  /*
    Read method to solve from  command line.  There seem to be several
    solvers that address a  problem of finding solvent distribution in
    external field  given the direct correlation function  of the pure
    solvent. This is a common entry  point for all of them. The actual
    algorith will be affected by --closure HNC/PY/KH and --snes-solver
    newton/picard/jager:
  */
  if (bgy3d_getopt_test ("--hnc"))
    {
      if (bgy3d_getopt_test ("--solute"))
        solver = HNC3d_solute_solve; /* solute/solvent by HNC */
      else
        solver = HNC3d_solvent_solve; /* pure solvent by HNC */
    }

  if (bgy3d_getopt_test ("--bgy"))
    {
      if (bgy3d_getopt_test ("--solute"))
        solver =  BGY3d_solute_solve; /* solvent/solute by BGY */
      else
        solver =  BGY3d_solvent_solve; /* pure solvent by BGY */
    }

  /* This one may only work for 3-site symmetric H2O-like solvents: */
  if (bgy3d_getopt_test ("--BGY3Site"))
    solver =  BGY3d_solvent_solve_h2o;

  if (solver)
    {
      local Vec g_ini;
      /* load initial configuration from file ??? */
      if (bgy3d_getopt_test ("--load"))
        {
          g_ini = bgy3d_vec_load ("g.bin");
          PetscPrintf (PETSC_COMM_WORLD, "g_ini loaded from file \"g.bin\".\n");
        }
      else
        g_ini = NULL;

      local Vec g = solver (&PD, g_ini);

      /* computation time measurement end point*/
      MPI_Barrier (PETSC_COMM_WORLD);
      mpi_stop = MPI_Wtime();
      PetscPrintf (PETSC_COMM_WORLD, "Total computation time: %.4f s\n",
                   mpi_stop - mpi_start);

      /* Output result */
      if (g != NULL)
        {
          bgy3d_vec_save_ascii ("vec.m", g);

          /* save g to binary file */
          if (bgy3d_getopt_test ("--save"))
            {
              bgy3d_vec_save ("g.bin", g);
              PetscPrintf (PETSC_COMM_WORLD, "Result written to file \"g.bin\".\n");
            }

          vec_destroy (&g);
        }
  }
  else
    PetscPrintf (PETSC_COMM_WORLD,
                 "Please choose one of: -BGY2site or -BGYM2site!\n");

  int ierr = PetscFinalize();CHKERRQ(ierr);

  /*
   * Original version of BGY3D  executable returned 1. Make interprets
   * non-zero codes as a failure, so this was changed to 0:
   */
  return 0;
}
#endif  /* else of ifdef WITH_GUILE */

#ifdef DEBUG
/* starts the gdb debugger for each process */
static int start_debugger (void)
{
  pid_t pid;
  char s[100];

  pid = getpid();

  sprintf(s, "xterm -e gdb /mount/kamino/jager/work/BGY3d/test/bgy3d %d &", pid);

  return system(s);
}
#endif
