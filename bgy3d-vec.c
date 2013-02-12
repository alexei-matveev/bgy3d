/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-vec.h"

Vec bgy3d_vec_create (const DA da)
{
  Vec g;
  DACreateGlobalVector (da, &g);
  return g;
}

void bgy3d_vec_destroy (Vec g)
{
  VecDestroy (g);
}

void bgy3d_vec_create1 (const DA da, int m, Vec g[m])
{
  for (int i = 0; i < m; i++)
    g[i] = bgy3d_vec_create (da);
}

void bgy3d_vec_destroy1 (int m, Vec g[m])
{
  for (int i = 0; i < m; i++)
    bgy3d_vec_destroy (g[i]);
}

/* Allocates g[m][m] with g[j][i] being aliased to g[i][j]: */
void bgy3d_vec_create2 (const DA da, int m, Vec g[m][m])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      g[j][i] = g[i][j] = bgy3d_vec_create (da);
}

void bgy3d_vec_destroy2 (int m, Vec g[m][m])
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      bgy3d_vec_destroy (g[i][j]);
}

static void shape (const DA da, int *NI, int *NJ, int *NK)
{
  int dim;
  DAGetInfo (da, &dim, NI, NJ, NK,
             NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  assert (dim == 3);
}


/*
  Does the mixing:

    cur := cur + a * (new - cur)

         = a * new + (1 - a) * cur

  Returns the norm of the difference |new - cur|.
 */
real bgy3d_vec_mix (Vec cur, Vec new, real a, Vec work)
{
  /* work := new - cur, in two steps: */
  VecCopy (new, work);
  VecAXPBY (work, -1.0, 1.0, cur);

  /* Norm of the change: */
  real norm;
  VecNorm (work, NORM_INFINITY, &norm);

  /* cur' = cur + a * (new - cur) */
  VecAXPBY (cur, a, 1.0, work);

  return norm;
}


/* This one is supposed to save enough meta-info (such as distribution
   pattern, dimensions) to recover the vector from scratch: */
void bgy3d_vec_save (const char file[], const Vec vec)
{
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_WRITE, &viewer);
  VecView (vec, viewer);
  PetscViewerDestroy (viewer);
}


/* This one  returns a newly  created vector: */
Vec bgy3d_vec_load (const char file[])
{
  Vec vec;                      /* new one */
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
  VecLoad (viewer, VECMPI, &vec); /* creates it */
  PetscViewerDestroy (viewer);

  return vec;
}


/* This one takes an allocated vector  and fills it with the data read
   from from disk. Pass a valid vector here: */
void bgy3d_vec_read (const char file[], Vec vec)
{
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
  VecLoadIntoVector (viewer, vec);
  PetscViewerDestroy (viewer);
}


/* This one will not save much of a meta-info: */
void bgy3d_vec_save_ascii (const char file[], const Vec vec)
{
  PetscViewer viewer;

  PetscViewerASCIIOpen (PETSC_COMM_WORLD, file, &viewer);
  /* PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_MATLAB); */
  /* PetscViewerSetFormat (viewer, PETSC_VIEWER_DEFAULT); */
  /* PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_VTK); */
  PetscViewerSetFormat (viewer, PETSC_VIEWER_ASCII_COMMON);
  VecView (vec, viewer);
  PetscViewerDestroy (viewer);
}


void bgy3d_vec_save_ascii1 (const char *format, int m, const Vec vec[m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Writing files...");
  for (int i = 0; i < m; i++)
    {
      char name[20 + strlen (format)];
      snprintf (name, sizeof name, format, i);
      bgy3d_vec_save_ascii (name, vec[i]);
    }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}


void bgy3d_vec_save_ascii2 (const char *format, int m, /* const */ Vec vec[m][m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Writing files...");
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        assert (vec[j][i] == vec[i][j]);

        char name[20 + strlen (format)];
        snprintf (name, sizeof name, format, j, i); /* ji as in g01.bin */

        bgy3d_vec_save_ascii (name, vec[i][j]);
      }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}


/*
  Read  solvent-solvent pair  distributions  g2[][] from  disk into  a
  pre-allocated   global   (distributed)   vectors.   Note   that   in
  bgy3d_vec_load(), the  type of  the Vec will  depend on  the on-disk
  data  which may  or  may not  be  compatible to  ones  used in  this
  run. That is why we use bgy3d_vec_read() here instead.

  Format is e.g. "g%d%d.bin"
*/
void bgy3d_vec_read2 (const char *format, int m, /* const */ Vec g2[m][m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Loading binary g2 files...");
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        assert (g2[j][i] == g2[i][j]);

        char name[20];
        snprintf (name, sizeof name, format, j, i); /* ji as in g01.bin */

        bgy3d_vec_read (name, g2[i][j]);
      }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}


/* Save  solvent-solvent  pair  distributions  g2[][]  to  disk.   See
   bgy3d_vec_read2(). */
void bgy3d_vec_save2 (const char *format, int m, /* const */ Vec g2[m][m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Writing binary g2 files...");
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        assert (g2[j][i] == g2[i][j]);

        char name[20];
        snprintf (name, sizeof name, format, j, i); /* ji as in g01.bin */

        bgy3d_vec_save (name, g2[i][j]);
      }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}

/* Format is e.g. "g%d.bin" */
void bgy3d_vec_save1 (const char *format, int m, const Vec g[m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Writing binary g1 files...");
  for (int i = 0; i < m; i++)
    {
      char name[20];
      snprintf (name, sizeof name, format, i);
      bgy3d_vec_save (name, g[i]);
    }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}


void bgy3d_vec_read1 (const char *format, int m, const Vec g[m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Loading binary files...");
  for (int i = 0; i < m; i++)
    {
      char name[20];
      snprintf (name, sizeof name, format, i);
      bgy3d_vec_read (name, g[i]);
    }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}

/* Fills Vec  g2 with  3D distribution derived  from the 1D  g(r) data
   from the disk.  Here g2 should be a valid allocated vector. */
static void vec_read_radial (const State *BHD, const char *filename, Vec g2)
{
  DA da;
  FILE *fp;
  real *xg, *g;
  real r[3], r_s, h[3], interval[2];
  int index=0;
  int x[3], n[3], i[3], k;
  PetscScalar ***g2_vec;

  da = BHD->da;
  const ProblemData *PD = BHD->PD;
  FOR_DIM
    h[dim] = PD->h[dim];

  interval[0] = PD->interval[0];
  interval[1] = PD->interval[1];

  /* read file */
  fp = fopen(filename, "r");
  if(fp==NULL)
    {
      PetscPrintf(PETSC_COMM_WORLD,"Could not open file %s.\n", filename);
      exit(1);
    }
  xg= (real*) malloc(sizeof(*xg));
  g= (real*) malloc(sizeof(*g));

  while( fscanf(fp,"%lf %lf", &xg[index], &g[index]) == 2)
    {
      index++;
      xg= (real*) realloc(xg, (index+1)*sizeof(*xg));
      g= (real*) realloc(g, (index+1)*sizeof(*g));
    }
  index--;
  PetscPrintf(PETSC_COMM_WORLD,"Read %d lines from file %s, x=[%f,%f]\n", index+1, filename,
              xg[0], xg[index]);
  fclose(fp);


  /* interpolate to 3d grid */

  /* Get local portion of the grid */
  DAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  DAVecGetArray(da, g2, &g2_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          FOR_DIM
            r[dim] = i[dim]*h[dim]+interval[0];
          r_s = sqrt( SQR(r[0])+SQR(r[1])+SQR(r[2]) );

          /* find x in array */
          for(k=0; k<=index; k++)
            {
            if(r_s<xg[k])
              break;
            }
          /* in area where r_s < xg[0] */
          if(k==0)
            g2_vec[i[2]][i[1]][i[0]] = 0;
          /* in area where r_s >= xg[index],
           * fill the vector with 1.0
           */
          else if(k>=index)
            g2_vec[i[2]][i[1]][i[0]] = 1.0;
          else
            /* Linear interpolation here:
             *                       x - x0
             * y = y0 + (y1 - y0) * --------
             *                      x1 - x0
             */
            g2_vec[i[2]][i[1]][i[0]] = g[k] + (r_s-xg[k])*(g[k+1]-g[k])/(xg[k+1]-xg[k]);

          /* FIXME: why do this?*/
          if(g2_vec[i[2]][i[1]][i[0]]<0)
            g2_vec[i[2]][i[1]][i[0]]=0;
        }
  DAVecRestoreArray(da, g2, &g2_vec);

  free(xg);
  free(g);
}


/*
  Was  used  for CS2  where  the g2  comes  from  MM simulations.   By
  convention bgy3d_vec_read* expects the storage for the vectors to be
  allocated.
*/
void bgy3d_vec_read_radial2 (const State *BHD,
                             const char *format, int m, /* const */ Vec g2[m][m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Loading radial g2 files...\n");
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        assert (g2[j][i] == g2[i][j]);

        char name[20];
        snprintf (name, sizeof name, format, j, i); /* ji as in g01.bin */

        vec_read_radial (BHD, name, g2[i][j]);
      }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}

void bgy3d_vec_moments (const DA da, Vec v,
                        real *q, real *x, real *y, real *z)
{
  /* Historically the grid origin is at 0.5 N[]: */
  int N[3];
  shape (da, &N[0], &N[1], &N[2]);

  real m0 (real v, int i, int j, int k)
  {
    (void) i; (void) j; (void) k;
    return v * 1;
  }
  real mx (real v, int i, int j, int k)
  {
    (void) j; (void) k;
    return v * (i - 0.5 * N[0]);
  }
  real my (real v, int i, int j, int k)
  {
    (void) i; (void) k;
    return v * (j - 0.5 * N[1]);
  }
  real mz (real v, int i, int j, int k)
  {
    (void) i; (void) j;
    return v * (k - 0.5 * N[2]);
  }

  *q = bgy3d_vec_integrate (da, m0, v);
  *x = bgy3d_vec_integrate (da, mx, v);
  *y = bgy3d_vec_integrate (da, my, v);
  *z = bgy3d_vec_integrate (da, mz, v);
}
