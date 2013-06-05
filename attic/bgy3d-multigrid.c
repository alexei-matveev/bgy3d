/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/*==========================================================*/
/*  $Id: bgy3d_MG.c,v 1.1 2007-07-31 17:12:33 jager Exp $ */
/*==========================================================*/

#include "bgy3d.h"
#include "bgy3d-multigrid.h"

#ifdef L_BOUNDARY_MG
/* Initialize M-Matrix with appropriate stencil */
PetscErrorCode  InitializeLaplaceMatrix(DMMG dmmg,Mat J,Mat M)
{
  State *BHD;
  PData PD;
  DA da;
  int x[3], n[3], i[3], N[3], dim, border;
  MatStencil col[3],row;
  PetscScalar v[3], vb=1.0;
  real h[3], L, zpad;

  PetscPrintf(PETSC_COMM_WORLD,"Assembling Matrix...");

  BHD = (State*) dmmg->user;

  da =  (DA)dmmg->dm;
  PD = BHD->PD;
  MatZeroEntries(M);
  L= PD->interval[1]-PD->interval[0];
  zpad = BHD->zpad;

  DAGetInfo(da,0,&(N[0]),&(N[1]),&(N[2]),0,0,0,0,0,0,0);
  PetscPrintf(PETSC_COMM_WORLD,"%d %d %d \n", N[0], N[1], N[2]);
  FOR_DIM
    h[dim] = L/N[dim];

  border = (int) ceil( ((PD->interval[1]-PD->interval[0])-(2.*zpad))/h[0]/2. );


  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));




  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  FOR_DIM
	    {
	      col[dim].i=i[0];
	      col[dim].j=i[1];
	      col[dim].k=i[2];
	      row.i=i[0];
	      row.j=i[1];
	      row.k=i[2];
	    }

	  /* Boundary */
	  if( i[0] <= border+1 || i[1] <= border+1 || i[2] <= border+1)
	    MatSetValuesStencil(M,1,&row,1,col+1,&vb,ADD_VALUES);
	  else if( i[0]>=N[0]-1-border || i[1]>=N[1]-1-border || i[2]>=N[2]-1-border)
	    MatSetValuesStencil(M,1,&row,1,col+1,&vb,ADD_VALUES);
	  else
	    {
	      FOR_DIM
		{
		  /* position in matrix */
		  switch(dim)
		    {
		    case 0: col[0].i -= 1; col[2].i += 1;break;
		    case 1: col[0].j -= 1; col[2].j += 1;break;
		    case 2: col[0].k -= 1; col[2].k += 1;break;
		    }
		  /* values to enter */
		  v[0]=1.0/SQR(h[dim]);
		  v[1]=-2.0/SQR(h[dim]);
		  v[2]=+1.0/SQR(h[dim]);


		  MatSetValuesStencil(M,1,&row,3,col,v,ADD_VALUES);
		  switch(dim)
		    {
		    case 0: col[0].i += 1; col[2].i -= 1;break;
		    case 1: col[0].j += 1; col[2].j -= 1;break;
		    case 2: col[0].k += 1; col[2].k -= 1;break;
		    }

		}
	    }
	}


  MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);

  PetscPrintf(PETSC_COMM_WORLD,"done.\n");
  return 0;
}

void CopyBoundary(State *BHD, Vec gfrom, Vec gto, real zpad)
{
  PData PD;
  DA da;
  int x[3], n[3], i[3], ic[3], dim, k, N[3], border;
  PetscScalar ***gfrom_vec, ***gto_vec;

  PD = BHD->PD;
  da = BHD->da;
  FOR_DIM
    N[dim] = PD->N[dim];

  VecSet(gto, 0.0);

  border = (int) ceil( ((PD->interval[1]-PD->interval[0])-(2.*zpad))/PD->h[0]/2. );

  /* Get local portion of the grid */
  DAGetCorners(da, &(x[0]), &(x[1]), &(x[2]), &(n[0]), &(n[1]), &(n[2]));

  DAVecGetArray(da, gfrom, &gfrom_vec);
  DAVecGetArray(da, gto, &gto_vec);


/*   FOR_DIM */
/*     { */
/*       ic[0]=dim; */
/*       ic[1]=(dim+1)%3; */
/*       ic[2]=(dim+2)%3; */
/*       i[ic[0]] = 0; */
/*       if( x[ic[0]] == 0 )  */
/* 	for( i[ic[1]]=x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++) */
/* 	  for( i[ic[2]]=x[ic[2]]; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++) */
/* 	    gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]]; */
/*       i[ic[0]] = N[ic[0]]-1; */
/*       if( x[ic[0]]+n[ic[0]] == N[ic[0]] )  */
/* 	for( i[ic[1]]=x[ic[1]]; i[ic[1]]<x[ic[1]]+n[ic[1]]; i[ic[1]]++) */
/* 	  for( i[ic[2]]=x[ic[2]]; i[ic[2]]<x[ic[2]]+n[ic[2]]; i[ic[2]]++) */
/* 	    gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]]; */
/*     } */

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
	{
	  if( ( i[0]<=border+1 || i[0]>=N[0]-1-border) ||
	      ( i[1]<=border+1 || i[1]>=N[1]-1-border) ||
	      ( i[2]<=border+1 || i[2]>=N[2]-1-border) )
	    gto_vec[i[2]][i[1]][i[0]] = gfrom_vec[i[2]][i[1]][i[0]];
	}


  DAVecRestoreArray(da, gfrom, &gfrom_vec);
  DAVecRestoreArray(da, gto, &gto_vec);
}

PetscErrorCode ComputeRHSDMMG(DMMG dmmg,Vec b)
{
  State *BHD;

  BHD = (State*) dmmg->user;

  VecCopy(BHD->pre, b);
  return 0;
}

void InitializeDMMGSolver(State *BHD)
{
  PC pc;
  DMMG *dmmg;
  int l;

  DMMGCreate(PETSC_COMM_WORLD,2,PETSC_NULL,&dmmg);
  DMMGSetDM(dmmg,(DM)BHD->da_dmmg);
  DADestroy(BHD->da_dmmg);
  for (l = 0; l < DMMGGetLevels(dmmg); l++)
    {
      DMMGSetUser(dmmg,l,(void*) BHD);
    }

  DMMGSetKSP(dmmg,ComputeRHSDMMG,InitializeLaplaceMatrix);


  BHD->dmmg = dmmg;
}


real ImposeLaplaceBoundary(State *BHD, Vec g, Vec b, Vec x, real zpad, int *iter)
{
  real mpi_start, mpi_stop;
  Vec sol;

  /* computation time measurement start point */
  MPI_Barrier( PETSC_COMM_WORLD);
  mpi_start = MPI_Wtime();

  /* Get boundary of g */
  CopyBoundary(BHD, g, BHD->pre, zpad);
  VecSet(x, 0.0);

  /* Solve Laplace */
  DMMGSolve(BHD->dmmg);

  sol = DMMGGetx(BHD->dmmg);
  VecCopy(sol, x);



  /* subtract solution from g */
  VecAXPY(g, -1.0, x);

  /* computation time measurement stop point */
  MPI_Barrier( PETSC_COMM_WORLD);
  mpi_stop = MPI_Wtime();

  return mpi_stop-mpi_start;

  /*  VecView(b,PETSC_VIEWER_STDERR_WORLD);   */
/*   exit(1);  */

}


#endif
