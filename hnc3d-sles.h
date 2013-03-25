/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */

#include <cblas.h>              /* cblas_dgemm/zgemm() */

/*
  FIXME:   is  there  a   C  interface   to  LAPACK?    Petsc  defines
  LAPACKgesv_()  macro  in  "petscblaslapack.h",  but it  resolves  to
  dgesv()  when PetscScalar  ==  real.  We  need  the complex  version
  though. The real version is here for completeness.
*/
void dgesv_ (const int *n, const int *nrhs,
             real *restrict a, const int *lda,
             int *ipiv,
             real *restrict b, const int *ldb,
             int *info);


void zgesv_ (const int *n, const int *nrhs,
             complex *restrict a, const int *lda,
             int *ipiv,
             complex *restrict b, const int *ldb,
             int *info);


/* Solve linear equations A X = B. As in LAPACK the matrix A is
   destroyed and the result is returned in B */
static inline void hnc3d_sles_dgesv (int m, real a[m][m], real b[m][m])
{
  int ipiv[m], info;

  /* B will be overwriten with the result, A will be overwritten
     with its factorization: */
  dgesv_ (&m, &m, (real*) a, &m, ipiv, (real*) b, &m, &info);
  assert (info == 0);
}


static inline void hnc3d_sles_zgesv (int m, complex a[m][m], complex b[m][m])
{
  int ipiv[m], info;

  /* B will be overwriten with the result, A will be overwritten
     with its factorization: */
  zgesv_ (&m, &m, (complex*) a, &m, ipiv, (complex*) b, &m, &info);
  assert (info == 0);
}


/* C := A * B. In fact, MM can be avoided. For small m it may not even
   pay off. */
static inline void hnc3d_sles_dgemm (int m,
                                     real a[m][m], /* in */
                                     real b[m][m], /* in */
                                     real c[m][m]) /* out */
{
  /*
    void cblas_dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                      const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                      const int K, const double alpha, const double *A,
                      const int lda, const double *B, const int ldb,
                      const double beta, double *C, const int ldc);
   */
  cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
               m, m, m, 1.0, (real*) a, m, (real*) b, m, 0.0, (real*) c, m);
}


static inline void hnc3d_sles_zgemm (int m,
                                     complex a[m][m], /* in */
                                     complex b[m][m], /* in */
                                     complex c[m][m]) /* out */
{
  /*
    void cblas_zgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                      const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                      const int K, const void *alpha, const void *A,
                      const int lda, const void *B, const int ldb,
                      const void *beta, void *C, const int ldc);

   */
  const complex one = 1.0, zero = 0.0;

  cblas_zgemm (CblasColMajor, CblasNoTrans, CblasNoTrans,
               m, m, m, &one, a, m, b, m, &zero, c, m);
}
