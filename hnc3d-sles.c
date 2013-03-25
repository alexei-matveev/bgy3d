/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include "bgy3d.h"
#include <complex.h>            /* before hnc3d-sles.h */
#include "hnc3d-sles.h"         /* also degsv_(), zgesv_() */


static void test_real ()
{
  /* both const: */
  real A[2][2] = {{3.0, 0.5},
                  {4.0, 1.0}};
  real B[2][2] = {{4.0, 8.0},
                  {3.0, 7.0}};

  const int m = 2;

  real a[m][m], x[m][m];

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        a[i][j] = A[i][j];
        x[i][j] = B[i][j];
      }

  hnc3d_sles_dgesv (m, a, x);

  real b[m][m];

  hnc3d_sles_dgemm (m, A, x, b);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        const real r = B[j][i] - b[j][i];
        printf ("diff = %f\n", r);
      }
}


static void test_complex ()
{
  /* both const: */
  complex A[2][2] = {{3.0 + I * 1.0, 0.5},
                     {4.0 + I * 1.0, 1.0}};
  complex B[2][2] = {{4.0, 8.0},
                     {3.0, 7.0 + I * 1.0}};

  const int m = 2;

  complex a[m][m], x[m][m];

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        a[i][j] = A[i][j];
        x[i][j] = B[i][j];
      }

  hnc3d_sles_zgesv (m, a, x);

  complex b[m][m];

  hnc3d_sles_zgemm (m, A, x, b);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      {
        const complex r = B[j][i] - b[j][i];
        printf ("diff = %f + %fi\n", creal (r), cimag (r));
      }
}


/* Not really used: */
void hnc3d_sles_test ()
{
  test_real ();
  test_complex ();
}
