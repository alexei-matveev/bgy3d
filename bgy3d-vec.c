/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include "bgy3d-vec.h"

#define assert_range(i, n) assert (likely (0 <= (i)) && likely ((i) < (n)));

/*
  Interpolate f(x) given the tabular values. The values in the x-table
  are

    xtab[i] = f(x ), for x = (i + 1/2) dx
                 i        i

  FIXME: precompute derivatives once things start becoming costly.
*/
static
real interp (real x, int n, const real xtab[n], real dx)
{
  const real ix = x / dx - 0.5; /* i(x), real! */

  /* Nearest table entry: */
  const int i = ix + 0.5;           /* rounding! */

  /* FIXME: Table should be complete enough: */
  assert_range (i, n);

  /* Here [ip, im] = [i + 1, i - 1] in most cases: */
  const int ip = (i < n - 2) ? i + 1 : i;
  const int im = (i > 0 + 1) ? i - 1 : i;
  assert_range (ip, n);         /* paranoya! */
  assert_range (im, n);         /* paranoya! */
  assert (ip - im > 0);         /* FIXME: fails for n < 2! */

  /* Estimate for the derivative df/di: */
  const real fi = (xtab[ip] - xtab[im]) / (ip - im);

  /* Note that  (ix - i) is  the *real* offset from  the nearest table
     entry: */
  return xtab[i] + fi * (ix - i);
}


void vec_rtab (const State *HD, int n, const real rtab[n], real dr,
                      Vec v) /* out */
{
  pure real f (real r)
  {
    return interp (r, n, rtab, dr);
  }
  vec_rmap (HD, f, v);
}


void vec_ktab (const State *HD, int n, const real ktab[n], real dk,
                      Vec v_fft) /* out */
{
  pure complex f (real k)
  {
    return interp (k, n, ktab, dk);
  }
  vec_kmap (HD, f, v_fft);
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
  PetscViewerDestroy (&viewer);
}


/* This one  returns a newly  created vector: */
Vec bgy3d_vec_load (const char file[])
{
  Vec vec;                      /* new one */
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
# if PETSC_VERSION < VERSION(3, 2)
  VecLoad (viewer, VECMPI, &vec); /* creates it */
#else
  /* FIXME: how to make it distributed? */
  VecCreate (PETSC_COMM_WORLD, &vec);
  VecLoad (vec, viewer);
#endif
  PetscViewerDestroy (&viewer);

  return vec;
}


/* This one takes an allocated vector  and fills it with the data read
   from from disk. Pass a valid vector here: */
void bgy3d_vec_read (const char file[], Vec vec)
{
  PetscViewer viewer;

  PetscViewerBinaryOpen (PETSC_COMM_WORLD, file, FILE_MODE_READ, &viewer);
  VecLoadIntoVector (viewer, vec);
  PetscViewerDestroy (&viewer);
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
  PetscViewerDestroy (&viewer);
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
   from the disk.  Here Vec g2 should be a valid allocated vector. */
void bgy3d_vec_read_radial (const DA da, const ProblemData *PD, const char *filename, Vec g2)
{
  FILE *fp;
  real *xg, *g;
  real r[3], r_s;
  int index=0;
  int x[3], n[3], i[3], k;
  PetscScalar ***g2_vec;

  const real *h = PD->h;        /* [3] */
  const real *L = PD->L;        /* [3] */

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
  DMDAGetCorners(da, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  DMDAVecGetArray(da, g2, &g2_vec);

  /* loop over local portion of grid */
  for(i[2]=x[2]; i[2]<x[2]+n[2]; i[2]++)
    for(i[1]=x[1]; i[1]<x[1]+n[1]; i[1]++)
      for(i[0]=x[0]; i[0]<x[0]+n[0]; i[0]++)
        {
          FOR_DIM
            r[dim] = i[dim] * h[dim] - L[dim] / 2;

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
  DMDAVecRestoreArray(da, g2, &g2_vec);

  free(xg);
  free(g);
}


/*
  Was  used  for CS2  where  the g2  comes  from  MM simulations.   By
  convention bgy3d_vec_read* expects the storage for the vectors to be
  allocated.
*/
void bgy3d_vec_read_radial2 (const DA da, const ProblemData *PD,
                             const char *format, int m, /* const */ Vec g2[m][m])
{
  PetscPrintf (PETSC_COMM_WORLD, "Loading radial g2 files...\n");
  for (int i = 0; i < m; i++)
    for (int j = 0; j <= i; j++)
      {
        assert (g2[j][i] == g2[i][j]);

        char name[20];
        snprintf (name, sizeof name, format, j, i); /* ji as in g01.bin */

        bgy3d_vec_read_radial (da, PD, name, g2[i][j]);
      }
  PetscPrintf (PETSC_COMM_WORLD, "done.\n");
}

void bgy3d_vec_moments1 (const DA da, Vec v, real *x, real *y, real *z)
{
  /* Historically the grid origin is at 0.5 N[]: */
  int N[3];
  da_shape (da, N);

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

  *x = vec_integrate (da, mx, v);
  *y = vec_integrate (da, my, v);
  *z = vec_integrate (da, mz, v);
}

/*
 * Second moments
 */
void bgy3d_vec_moments2 (const DA da, Vec v,
                         real *xy, real *yz, real *zx,
                         real *z2, real *x2y2, real *r2)
{
  /* Historically the grid origin is at 0.5 N[]: */
  int N[3];
  da_shape (da, N);

  /* < x * y > */
  real mxy (real v, int i, int j, int k)
  {
    (void) k;
    return v * (i - 0.5 * N[0]) * (j - 0.5 * N[1]);
  }
  /* < y * z > */
  real myz (real v, int i, int j, int k)
  {
    (void) i;
    return v * (j - 0.5 * N[1]) * (k - 0.5 * N[2]);
  }
  /* < z * x > */
  real mzx (real v, int i, int j, int k)
  {
    (void) j;
    return v * (k - 0.5 * N[2]) * (i - 0.5 * N[0]);
  }
  /* < z^2 - 1 / 3 * r^2 > */
  real mz2 (real v, int i, int j, int k)
  {
    return v * (2 * SQR (k - 0.5 * N[2]) -
                SQR (i - 0.5 * N[0]) -
                SQR (j - 0.5 * N[1])) / 3.0;
  }
  /* < x^2 - y^2 > */
  real mx2y2 (real v, int i, int j, int k)
  {
    (void) k;
    return v * (SQR (i - 0.5 * N[0]) - SQR (j - 0.5 * N[1]));
  }
  /* < x^2 + y^2 + z^2 > */
  real mr2 (real v, int i, int j, int k)
  {
    return v * (SQR (i - 0.5 * N[0]) +
                SQR (j - 0.5 * N[1]) +
                SQR (k - 0.5 * N[2]));
  }

  *xy   = vec_integrate (da, mxy, v);
  *yz   = vec_integrate (da, myz, v);
  *zx   = vec_integrate (da, mzx, v);
  *z2   = vec_integrate (da, mz2, v);
  *x2y2 = vec_integrate (da, mx2y2, v);
  *r2   = vec_integrate (da, mr2, v);
}


/*
  By scaling k-components of the FFT Vec v along each dimension by

    exp [i * (2πk/L) * L/2] = cos (πk) = (-1)^k

  we effectively  translate the center of  the grid to  the corner and
  vice  versa.   FIXME:  again,  since the  complex  array  descriptor
  contians different dimensions  we have to pass the  real grid shape,
  N[3], separately. Complex Vec v is intent(inout) here:
*/
void bgy3d_vec_fft_trans (const DA dc, const int N[static 3], Vec v)
{
  int x[3], n[3], i[3];

  /* Get local portion of the grid */
  DMDAGetCorners (dc, &x[0], &x[1], &x[2], &n[0], &n[1], &n[2]);

  /* Loop over local portion of grid: */
  complex ***v_;
  DMDAVecGetArray (dc, v, &v_);

  for (i[2] = x[2]; i[2] < x[2] + n[2]; i[2]++)
    for (i[1] = x[1]; i[1] < x[1] + n[1]; i[1]++)
      for (i[0] = x[0]; i[0] < x[0] + n[0]; i[0]++)
        {
          int ic[3];

          /* Take negative frequencies for i > N/2: */
          FOR_DIM
            ic[dim] = KFREQ (i[dim], N[dim]);

          /* phase shift factor for x=x+L/2 */
          const int sign = COSSIGN(ic[0]) * COSSIGN(ic[1]) * COSSIGN(ic[2]);

          v_[i[2]][i[1]][i[0]] *= sign;
        }
  DMDAVecRestoreArray (dc, v, &v_);
}
