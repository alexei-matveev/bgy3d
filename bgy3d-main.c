/*==========================================================*/
/*  $Id: bgy3d.c,v 1.48 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/


/*==========================================*/
/* Hinweise zur Faltung: */
/* g(x) muss immer um 0 herum bekannt sein. Also sollten auch */
/* die x_M um 0 herum liegen. */
/*===============================================*/

#include "bgy3d.h"
#include "bgy3d-getopt.h"
#include "bgy3dH2OS.h"          /* BGY3d_solve_H2O_2site */
#include "bgy3dH2O.h"           /* BGY3d_solve_2site */
#include "bgy3dmolecule.h"      /* BGY3d_solve_DiatomicAB */
#include "bgy3dH2ONewton.h"
#include "bgy3dH2OSNewton.h"
#include "hnc3d.h"
#include "bgy3ddiv.h"           /* BGY3dDiv_solve, BGY3dDiv_solve2 */
#include "bgy3dtest.h"          /* BGY3dDiv_test */
#include "bgy3dfourier.h"       /* BGY3dDiv_solve_Fourier */
#include "bgy3d-simple.h"       /* BGY3d_solve */

static void PData_CreateParallel (ProblemData *PD);
static int start_debugger (void);

static char helptext[] = "Solving BGY3d equation.\n";

int verbosity = 0;

typedef Vec Solver (ProblemData *PD, Vec g_ini);

int main (int argc, char **argv)
{
  ProblemData PD;
  int ierr, N=0;
  real beta=0.6061, rho=0.3, h=0.5, interval[2]={-5.0,5.0};
  real mpi_start, mpi_stop;
  int np;
  Vec g, g_ini;
  Solver *solver = NULL;

  verbosity=0;

  PetscInitialize( &argc, &argv, (char*)0, helptext);


  MPI_Comm_size(PETSC_COMM_WORLD, &np);
  PetscPrintf(PETSC_COMM_WORLD, "NP= %d\n", np);


  /* computation time measurement start point */
  MPI_Barrier( PETSC_COMM_WORLD);
  mpi_start = MPI_Wtime();

  /* Read the command  line options. Petsc insists on  keys having the
     leading dash, so keep them for the moment. */

  /* Grid points in 1 dimension */
  bgy3d_getopt_int ("-N", &N);

  /* inverse temperature */
  bgy3d_getopt_real ("-beta", &beta);

  /* Density */
  bgy3d_getopt_real ("-rho", &rho);

  /* set verbosity */
  bgy3d_getopt_int ("-v", &verbosity);

  /* mesh width */
  bgy3d_getopt_real ("-mesh", &h);

  /*=================================*/
  /* set Problem Data */
  /*=================================*/

  if(N==0)
    N= (int) ceil((interval[1]-interval[0])/h);
  else
    h = (interval[1]-interval[0])/N;
  FOR_DIM
    PD.N[dim] = N;
  PD.N3 = PD.N[0]*PD.N[1]*PD.N[2];
  FOR_DIM
    PD.h[dim] = h;
  PD.interval[0] = interval[0];
  PD.interval[1] = interval[1];
  PD.beta = beta;
  PD.rho  = rho;

  PData_CreateParallel(&PD);
  PetscViewerSetFormat(PETSC_VIEWER_STDERR_WORLD,PETSC_VIEWER_ASCII_MATLAB);
  PetscViewerSetFormat(PETSC_VIEWER_STDERR_SELF,PETSC_VIEWER_ASCII_MATLAB);
  /*==================================*/

  PetscPrintf(PETSC_COMM_WORLD, "Grid size N=%d %d %d\n",PD.N[0], PD.N[1], PD.N[2]);
  PetscPrintf(PETSC_COMM_WORLD, "Total dof N^3=%d\n",PD.N3);
  PetscPrintf(PETSC_COMM_WORLD, "Domain [%f %f]^3\n", interval[0], interval[1]);
  PetscPrintf(PETSC_COMM_WORLD, "h = %f\n", h);
  PetscPrintf(PETSC_COMM_WORLD, "beta = %f\n", beta);
  PetscPrintf(PETSC_COMM_WORLD, "rho = %f\n", rho);

  //PetscPrintf(PETSC_COMM_WORLD, "\tATTENTION: Factor 2 is included!!! But why???\n");

  //if(PD.id==1)
/*   start_debugger(); */
/*   sleep(5); */


  /* Read method to solve from command line */
  if (bgy3d_getopt_test ("-simple"))
    solver = BGY3d_solve;

  if (bgy3d_getopt_test ("-HNC"))
    solver = HNC3d_Solve_h;

  if (bgy3d_getopt_test ("-HNCNewton"))
    solver = HNC3dNewton2_solve;

  if (bgy3d_getopt_test ("-DIV"))
    solver = BGY3dDiv_solve2;

  if (bgy3d_getopt_test ("-DIV2"))
    solver = BGY3dDiv_solve;

  if (bgy3d_getopt_test ("-BGYTEST"))
    solver = BGY3dDiv_test;

  if (bgy3d_getopt_test ("-BGYFOURIER"))
    solver =  BGY3dDiv_solve_Fourier;

  if (bgy3d_getopt_test ("-BGYFOURIERTEST"))
    solver =  BGY3dDiv_solve_FourierTest;

  if (bgy3d_getopt_test ("-BGYCONVOLUTIONTEST"))
    solver =  BGY3d_Convolution_Test;

  if (bgy3d_getopt_test ("-BGYDIATOMIC"))
    solver =  BGY3d_solve_DiatomicAB;

  if (bgy3d_getopt_test ("-BGY2Site"))
    solver =  BGY3d_solve_2site;

  if (bgy3d_getopt_test ("-BGY3Site"))
    solver =  BGY3d_solve_3site;

  if (bgy3d_getopt_test ("-BGYM2Site"))
    solver =  BGY3dM_solve_H2O_2site;

  if (bgy3d_getopt_test ("-BGYM3Site"))
    solver =  BGY3dM_solve_H2O_3site;

  if (bgy3d_getopt_test ("-BGYH2ONEWTON"))
    solver = BGY3d_SolveNewton_H2O;

  if (bgy3d_getopt_test ("-BGYH2OSNEWTON"))
    solver = BGY3d_SolveNewton_H2OS;

  if (bgy3d_getopt_test ("-BGYH2OSFNEWTON"))
    solver = BGY3d_SolveNewton_H2OSF;

  if(solver) {
      /* load initial configuration from file ??? */
      if (bgy3d_getopt_test ("-load")) {
          bgy3d_load_vec ("g.bin", &g_ini);
          PetscPrintf(PETSC_COMM_WORLD,"g_ini loaded from file \"g.bin\".\n");
      }
      else
          g_ini = PETSC_NULL;

      g = solver (&PD, g_ini);

      /* computation time measurement end point*/
      MPI_Barrier( PETSC_COMM_WORLD);
      mpi_stop = MPI_Wtime();
      PetscPrintf(PETSC_COMM_WORLD,"Total computation time: %.4f s\n",
                  mpi_stop-mpi_start);



      /* Output result */
      if( g != PETSC_NULL) {
          bgy3d_save_vec_ascii ("vec.m", g);

          /* save g to binary file */
          if (bgy3d_getopt_test ("-save")) {
              bgy3d_save_vec ("g.bin", g);
	      PetscPrintf(PETSC_COMM_WORLD,"Result written to file \"g.bin\".\n");
          }

	  VecDestroy(g);
      }

  }
  else
    PetscPrintf(PETSC_COMM_WORLD, "Please choose one of: -BGY2site or -BGYM2site!\n");

  ierr = PetscFinalize();CHKERRQ(ierr);

  /*
   * Original version of BGY3D  executable returned 1. Make interprets
   * non-zero codes as a failure, so this was changed to 0:
   */
  return 0;
}

static void PData_CreateParallel (ProblemData *PD)
{

  MPI_Comm_size(PETSC_COMM_WORLD, &(PD->np));
  MPI_Comm_rank(PETSC_COMM_WORLD, &(PD->id));
}

/* starts the gdb debugger for each process */
static int start_debugger (void)
{
  pid_t pid;
  char s[100];

  pid = getpid();

  sprintf(s, "xterm -e gdb /mount/kamino/jager/work/BGY3d/test/bgy3d %d &", pid);

  return system(s);
}
