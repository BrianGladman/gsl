/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Implementation of gsl_blas_raw interface which
 * defers to a CBLAS conformant interface.
 */
#include "gsl_blas_raw.h"

/*
 * ===========================================================================
 * level 1 BLAS functions
 * ===========================================================================
 */

float  gsl_blas_raw_sdsdot (size_t N,
                            float alpha,
                            const float X[], size_t incX,
                            const float Y[], size_t incY)
{
  return cblas_sdsdot(N, alpha, X, incX, Y, incY);
}


double gsl_blas_raw_dsdot (size_t N,
                           const float X[], size_t incX,
                           const float Y[], size_t incY)
{
  return cblas_dsdot(N, X, incX, Y, incY);
}


float  gsl_blas_raw_sdot (size_t N,
                          const float X[], size_t incX,
                          const float Y[], size_t incY)
{
  return cblas_sdot(N, X, incX, Y, incY);
}

                        
double gsl_blas_raw_ddot (size_t N,
                          const double X[], size_t incX,
                          const double Y[], size_t incY)
{
  return cblas_ddot(N, X, incX, Y, incY);
}


void gsl_blas_raw_cdotu (size_t N,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_float dotu)
{
  cblas_cdotu_sub(N, X, incX, Y, incY, dotu);
}


void gsl_blas_raw_cdotc (size_t N,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_float dotc)
{
  cblas_cdotc_sub(N, X, incX, Y, incY, dotc);
}


void gsl_blas_raw_zdotu (size_t N,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
                         gsl_complex_packed dotu)
{
  cblas_zdotu_sub(N, X, incX, Y, incY, dotu);
}


void gsl_blas_raw_zdotc (size_t N,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
                         gsl_complex_packed dotc)
{
  cblas_zdotc_sub(N, X, incX, Y, incY, dotc);
}


float  gsl_blas_raw_snrm2  (size_t N, const float  X[], size_t incX)
{
  return cblas_snrm2(N, X, incX);
}


double gsl_blas_raw_dnrm2  (size_t N, const double X[], size_t incX)
{
  return cblas_dnrm2(N, X, incX);
}


float  gsl_blas_raw_scnrm2 (size_t N, const gsl_const_complex_packed_array_float X, size_t incX)
{
  return cblas_scnrm2(N, X, incX);
}


double gsl_blas_raw_dznrm2 (size_t N, const gsl_const_complex_packed_array X, size_t incX)
{
  return cblas_dznrm2(N, X, incX);
}


float  gsl_blas_raw_sasum  (size_t N, const float X[], size_t incX)
{
  return cblas_sasum(N, X, incX);
}


double gsl_blas_raw_dasum  (size_t N, const double X[], size_t incX)
{
  return cblas_dasum(N, X, incX);
}


float  gsl_blas_raw_scasum (size_t N, const gsl_const_complex_packed_array_float X, size_t incX)
{
  return cblas_scasum(N, X, incX);
}


double gsl_blas_raw_dzasum (size_t N, const gsl_const_complex_packed_array X, size_t incX)
{
  return cblas_dzasum(N, X, incX);
}


CBLAS_INDEX_t gsl_blas_raw_isamax (size_t N, const float X[], size_t incX)
{
  return cblas_isamax(N, X, incX);
}


CBLAS_INDEX_t gsl_blas_raw_idamax (size_t N, const double X[], size_t incX)
{
  return cblas_idamax(N, X, incX);
}


CBLAS_INDEX_t gsl_blas_raw_icamax (size_t N, const gsl_const_complex_packed_array_float X, size_t incX)
{
  return cblas_icamax(N, X, incX);
}


CBLAS_INDEX_t gsl_blas_raw_izamax (size_t N, const gsl_const_complex_packed_array X, size_t incX)
{
  return cblas_izamax(N, X, incX);
}


void gsl_blas_raw_sswap (size_t N,
                         float X[], size_t incX,
                         float Y[], size_t incY)
{
  cblas_sswap(N, X, incX, Y, incY);
}


void gsl_blas_raw_dswap (size_t N,
                         double X[], size_t incX,
                         double Y[], size_t incY)
{
  cblas_dswap(N, X, incX, Y, incY);
}


void gsl_blas_raw_cswap (size_t N,
                         gsl_complex_packed_array_float X, size_t incX,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_cswap(N, X, incX, Y, incY);
}


void gsl_blas_raw_zswap (size_t N,
                         gsl_complex_packed_array X, size_t incX,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zswap(N, X, incX, Y, incY);
}


void gsl_blas_raw_scopy (size_t N,
                         const float X[], size_t incX,
                         float Y[], size_t incY)
{
  cblas_scopy(N, X, incX, Y, incY);
}


void gsl_blas_raw_dcopy (size_t N,
                         const double X[], size_t incX,
                         double Y[], size_t incY)
{
  cblas_dcopy(N, X, incX, Y, incY);
}


void gsl_blas_raw_ccopy (size_t N,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_ccopy(N, X, incX, Y, incY);
}


void gsl_blas_raw_zcopy (size_t N,
                         const gsl_const_complex_packed_array X, size_t incX,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zcopy(N, X, incX, Y, incY);
}


void gsl_blas_raw_saxpy (size_t N,
                         float alpha,
                         const float X[], size_t incX,
                         float Y[], size_t incY)
{
  cblas_saxpy(N, alpha, X, incX, Y, incY);
}


void gsl_blas_raw_daxpy (size_t N,
                         double alpha,
                         const double X[], size_t incX, 
                         double Y[], size_t incY)
{
  cblas_daxpy(N, alpha, X, incX, Y, incY);
}


void gsl_blas_raw_caxpy (size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_caxpy(N, alpha, X, incX, Y, incY);
}


void gsl_blas_raw_zaxpy (size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array X, size_t incX,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zaxpy(N, alpha, X, incX, Y, incY);
}


void gsl_blas_raw_srotg (float a[], float b[], float c[], float s[])
{
  cblas_srotg(a, b, c, s);
}


void gsl_blas_raw_drotg (double a[], double b[], double c[], double s[])
{
  cblas_drotg(a, b, c, s);
}


void gsl_blas_raw_srotmg (float d1[], float d2[], float b1[],
                          float b2, float P[])
{
  cblas_srotmg(d1, d2, b1, b2, P);
}


void gsl_blas_raw_drotmg (double d1[], double d2[], double b1[],
                          double b2, double P[])
{
  cblas_drotmg(d1, d2, b1, b2, P);
}


void gsl_blas_raw_srot (size_t N,
                        float X[], size_t incX,
                        float Y[], size_t incY,
                        float c, float s)
{
  cblas_srot(N, X, incX, Y, incY, c, s);
}


void gsl_blas_raw_drot (size_t N,
                        double X[], size_t incX,
                        double Y[], size_t incY,
                        const double c, const double s)
{
  cblas_drot(N, X, incX, Y, incY, c, s);
}


void gsl_blas_raw_srotm (size_t N,
                         float X[], size_t incX,
                         float Y[], size_t incY,
                         const float P[])
{
  cblas_srotm(N, X, incX, Y, incY, P);
}


void gsl_blas_raw_drotm (size_t N,
                         double X[], size_t incX,
                         double Y[], size_t incY,
                         const double P[])
{
  cblas_drotm(N, X, incX, Y, incY, P);
}


void gsl_blas_raw_sscal  (size_t N, float  alpha, float  X[], size_t incX)
{
  cblas_sscal(N, alpha, X, incX);
}


void gsl_blas_raw_dscal  (size_t N, double alpha, double X[], size_t incX)
{
  cblas_dscal(N, alpha, X, incX);
}


void gsl_blas_raw_cscal  (size_t N, const gsl_const_complex_packed_float alpha, gsl_complex_packed_array_float X, size_t incX)
{
  cblas_cscal(N, alpha, X, incX);
}


void gsl_blas_raw_zscal  (size_t N, const gsl_const_complex_packed alpha, gsl_complex_packed_array X, size_t incX)
{
  cblas_zscal(N, alpha, X, incX);
}


void gsl_blas_raw_csscal (size_t N, float  alpha, gsl_complex_packed_array_float X, size_t incX)
{
  cblas_csscal(N, alpha, X, incX);
}


void gsl_blas_raw_zdscal (size_t N, double alpha, gsl_complex_packed_array X, size_t incX)
{
  cblas_zdscal(N, alpha, X, incX);
}


/*
 * ===========================================================================
 * level 2 BLAS
 * ===========================================================================
 */

 
/* GEMV */

void gsl_blas_raw_sgemv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N,
                         float alpha,
                         const float A[], int lda,
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
  cblas_sgemv(CblasRowMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dgemv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
  cblas_dgemv(CblasRowMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_cgemv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_cgemv(CblasRowMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zgemv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zgemv(CblasRowMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


/* GBMV */

void gsl_blas_raw_sgbmv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         float alpha,
                         const float A[], int lda,
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
  cblas_sgbmv(CblasRowMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dgbmv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
  cblas_dgbmv(CblasRowMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_cgbmv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_cgbmv(CblasRowMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zgbmv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zgbmv(CblasRowMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


/* TRMV */

void gsl_blas_raw_strmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const float A[], int lda,
                         float X[], size_t incX)
{
  cblas_strmv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_dtrmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const double A[], int lda,
                         double X[], size_t incX)
{
  cblas_dtrmv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ctrmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX)
{
  cblas_ctrmv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ztrmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX)
{
  cblas_ztrmv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


/* TBMV */

void gsl_blas_raw_stbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const float A[], int lda,
                         float X[], size_t incX)
{
  cblas_stbmv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_dtbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const double A[], int lda,
                         double X[], size_t incX)
{
  cblas_dtbmv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ctbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX)
{
  cblas_ctbmv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ztbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX)
{
  cblas_ztbmv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


/* TPMV */

void gsl_blas_raw_stpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const float Ap[],
                         float X[], size_t incX)
{
  cblas_stpmv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_dtpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const double Ap[],
                         double X[], size_t incX)
{
  cblas_dtpmv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ctpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array_float Ap,
                         gsl_complex_packed_array_float X, size_t incX)
{
  cblas_ctpmv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ztpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array Ap,
                         gsl_complex_packed_array X, size_t incX)
{
  cblas_ztpmv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


/* TRSV */

void gsl_blas_raw_strsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const float A[], int lda,
                         float X[], size_t incX)
{
  cblas_strsv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_dtrsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const double A[], int lda,
                         double X[], size_t incX)
{
  cblas_dtrsv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ctrsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX)
{
  cblas_ctrsv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ztrsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX)
{
  cblas_ztrsv(CblasRowMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


/* TBSV */

void gsl_blas_raw_stbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const float A[], int lda,
                         float X[], size_t incX)
{
  cblas_stbsv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_dtbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const double A[], int lda,
                         double X[], size_t incX)
{
  cblas_dtbsv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ctbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX)
{
  cblas_ctbsv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ztbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
                         const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX)
{
  cblas_ztbsv(CblasRowMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


/* TPSV */

void gsl_blas_raw_stpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const float Ap[],
                         float X[], size_t incX)
{
  cblas_stpsv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_dtpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const double Ap[],
                         double X[], size_t incX)
{
  cblas_dtpsv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ctpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array_float Ap,
                         gsl_complex_packed_array_float X, size_t incX)
{
  cblas_ctpsv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ztpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const gsl_const_complex_packed_array Ap,
                         gsl_complex_packed_array X, size_t incX)
{
  cblas_ztpsv(CblasRowMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


/* SYMV */

void gsl_blas_raw_ssymv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         float alpha,
                         const float A[], int lda,
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
  cblas_ssymv(CblasRowMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dsymv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
  cblas_dsymv(CblasRowMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


/* SBMV */

void gsl_blas_raw_ssbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
                         float alpha,
                         const float A[], int lda, 
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
  cblas_ssbmv(CblasRowMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dsbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
  cblas_dsbmv(CblasRowMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


/* SPMV */

void gsl_blas_raw_sspmv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         float alpha,
                         const float Ap[],
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
  cblas_sspmv(CblasRowMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dspmv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         double alpha,
                         const double Ap[],
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
  cblas_dspmv(CblasRowMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


/* GER */

void gsl_blas_raw_sger (size_t M, size_t N,
                        float alpha,
                        const float X[], size_t incX,
                        const float Y[], size_t incY,
                        float A[], int lda)
{
  cblas_sger(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_dger (size_t M, size_t N,
                        double alpha,
                        const double X[], size_t incX,
                        const double Y[], size_t incY,
                        double A[], int lda)
{
  cblas_dger(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


/* SYR */

void gsl_blas_raw_ssyr (CBLAS_UPLO_t Uplo,
                        size_t N,
                        float alpha,
                        const float X[], size_t incX,
                        float A[], int lda)
{
  cblas_ssyr(CblasRowMajor, Uplo, N,alpha, X, incX, A, lda);
}


void gsl_blas_raw_dsyr (CBLAS_UPLO_t Uplo,
                        size_t N,
                        double alpha,
                        const double X[], size_t incX,
                        double A[], int lda)
{
  cblas_dsyr(CblasRowMajor, Uplo, N, alpha, X, incX, A, lda);
}


/* SPR */

void gsl_blas_raw_sspr (CBLAS_UPLO_t Uplo,
                        size_t N,
                        float alpha,
                        const float X[], size_t incX,
                        float Ap[])
{
  cblas_sspr(CblasRowMajor, Uplo, N, alpha, X, incX, Ap);
}


void gsl_blas_raw_dspr (CBLAS_UPLO_t Uplo,
                        size_t N,
                        double alpha,
                        const double X[], size_t incX,
                        double Ap[])
{
  cblas_dspr(CblasRowMajor, Uplo, N, alpha, X, incX, Ap);
}


/* SYR2 */

void gsl_blas_raw_ssyr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         float alpha,
                         const float X[], size_t incX, 
                         const float Y[], size_t incY,
                         float A[], int lda)
{
  cblas_ssyr2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_dsyr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         double alpha,
                         const double X[], size_t incX,
                         const double Y[], size_t incY,
                         double A[], int lda)
{
  cblas_dsyr2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


/* SPR2 */

void gsl_blas_raw_sspr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         float alpha,
                         const float X[], size_t incX,
                         const float Y[], size_t incY,
                         float A[])
{
  cblas_sspr2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, A);
}


void gsl_blas_raw_dspr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         double alpha,
                         const double X[], size_t incX,
                         const double Y[], size_t incY,
                         double A[])
{
  cblas_sspr2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, A);
}


/* HEMV */

void gsl_blas_raw_chemv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_chemv(CblasRowMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zhemv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zhemv(CblasRowMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


/* HBMV */

void gsl_blas_raw_chbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_chbmv(CblasRowMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zhbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zhbmv(CblasRowMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


/* HPMV */

void gsl_blas_raw_chpmv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float Ap,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  cblas_chpmv(CblasRowMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zhpmv (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array Ap,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
  cblas_zhpmv(CblasRowMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


/* GERU */

void gsl_blas_raw_cgeru (size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float A, int lda)
{
  cblas_cgeru(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_zgeru (size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array A, int lda)
{
  cblas_zgeru(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


/* GERC */

void gsl_blas_raw_cgerc (size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float A, int lda)
{
  cblas_cgerc(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_zgerc (size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array A, int lda)
{
  cblas_zgerc(CblasRowMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


/* HER */

void gsl_blas_raw_cher (CBLAS_UPLO_t Uplo,
                        size_t N,
                        float alpha,
                        const gsl_const_complex_packed_array_float X, size_t incX,
                        gsl_complex_packed_array_float A, int lda)
{
  cblas_cher(CblasRowMajor, Uplo, N, X, incX, A, lda);
}


void gsl_blas_raw_zher (CBLAS_UPLO_t Uplo,
                        size_t N,
                        double alpha,
                        const gsl_const_complex_packed_array X, size_t incX,
                        gsl_complex_packed_array A, int lda)
{
  cblas_zher(CblasRowMajor, Uplo, N, X, incX, A, lda);
}


/* HPR */

void gsl_blas_raw_chpr (CBLAS_UPLO_t Uplo,
                        size_t N,
                        float alpha,
                        const gsl_const_complex_packed_array_float X, size_t incX,
                        gsl_complex_packed_array_float A)
{
  cblas_chpr(CblasRowMajor, Uplo, N, alpha, X, incX, A);
}


void gsl_blas_raw_zhpr (CBLAS_UPLO_t Uplo,
                        size_t N,
                        double alpha,
                        const gsl_const_complex_packed_array X, size_t incX,
                        gsl_complex_packed_array A)
{
  cblas_zhpr(CblasRowMajor, Uplo, N, alpha, X, incX, A);
}


/* HER2 */

void gsl_blas_raw_cher2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float A, int lda)
{
  cblas_cher2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_zher2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array A, int lda)
{
  cblas_zher2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


/* HPR2 */

void gsl_blas_raw_chpr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float Ap)
{
  cblas_chpr2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, Ap);
}


void gsl_blas_raw_zhpr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array Ap)
{
  cblas_zhpr2(CblasRowMajor, Uplo, N, alpha, X, incX, Y, incY, Ap);
}



/*
 * ===========================================================================
 * level 3 BLAS
 * ===========================================================================
 */


/* GEMM */

void gsl_blas_raw_sgemm (CBLAS_TRANSPOSE_t TransA,
                         CBLAS_TRANSPOSE_t TransB,
                         size_t M, size_t N, size_t K,
                         float alpha,
                         const float A[], int lda,
                         const float B[], int ldb,
                         float beta,
                         float C[], int ldc)
{
  cblas_sgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_dgemm (CBLAS_TRANSPOSE_t TransA,
                         CBLAS_TRANSPOSE_t TransB,
                         size_t M, size_t N, size_t K,
                         double alpha,
                         const double A[], int lda,
                         const double B[], int ldb,
                         double beta,
                         double C[], int ldc)
{
  cblas_dgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_cgemm (CBLAS_TRANSPOSE_t TransA,
                         CBLAS_TRANSPOSE_t TransB,
                         size_t M, size_t N, size_t K,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda, 
                         const gsl_const_complex_packed_array_float B, int ldb,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
  cblas_cgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zgemm (CBLAS_TRANSPOSE_t TransA,
                         CBLAS_TRANSPOSE_t TransB,
                         size_t M, size_t N, size_t K,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array B, int ldb,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array C, int ldc)
{
  cblas_zgemm(CblasRowMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* SYMM */

void gsl_blas_raw_ssymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo,
                         size_t M, size_t N,
                         float alpha,
                         const float A[], int lda,
                         const float B[], int ldb,
                         float beta,
                         float C[], int ldc)
{
  cblas_ssymm(CblasRowMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_dsymm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo,
                         size_t M, size_t N,
                         double alpha,
                         const double A[], int lda,
                         const double B[], int ldb,
                         double beta,
                         double C[], int ldc)
{
  cblas_dsymm(CblasRowMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);

}


void gsl_blas_raw_csymm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo,
                         size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_array_float B, int ldb,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
  cblas_csymm(CblasRowMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zsymm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo,
                         size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array B, int ldb,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array C, int ldc)
{
  cblas_zsymm(CblasRowMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* SYRK */

void gsl_blas_raw_ssyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                         size_t N, size_t K,
                         float alpha,
                         const float A[], int lda,
                         float beta,
                         float C[], int ldc)
{
  cblas_ssyrk(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_dsyrk (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t Trans,
                         size_t N, size_t K,
                         double alpha,
                         const double A[], int lda,
                         double beta,
                         double C[], int ldc)
{
  cblas_dsyrk(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_csyrk (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t Trans,
                         size_t N, size_t K,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
  cblas_csyrk(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_zsyrk (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t Trans,
                         size_t N, size_t K,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array C, int ldc)
{
  cblas_zsyrk(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


/* SYR2K */

void gsl_blas_raw_ssyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                          size_t N, size_t K,
                          float alpha,
                          const float A[], int lda,
                          const float B[], int ldb,
                          float beta,
                          float C[], int ldc)
{
  cblas_ssyr2k(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_dsyr2k (CBLAS_UPLO_t Uplo,
                          CBLAS_TRANSPOSE_t Trans,
                          size_t N, size_t K,
                          double alpha,
                          const double A[], int lda,
                          const double B[], int ldb,
                          double beta,
                          double C[], int ldc)
{
  cblas_dsyr2k(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_csyr2k (CBLAS_UPLO_t Uplo,
                          CBLAS_TRANSPOSE_t Trans,
                          size_t N, size_t K,
                          const gsl_const_complex_packed_float alpha,
                          const gsl_const_complex_packed_array_float A, int lda,
                          const gsl_const_complex_packed_array_float B, int ldb,
                          const gsl_const_complex_packed_float beta,
                          gsl_complex_packed_array_float C, int ldc)
{
  cblas_csyr2k(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zsyr2k (CBLAS_UPLO_t Uplo,
                          CBLAS_TRANSPOSE_t Trans,
                          size_t N, size_t K,
                          const gsl_const_complex_packed alpha,
                          const gsl_const_complex_packed_array A, int lda,
                          const gsl_const_complex_packed_array B, int ldb,
                          const gsl_const_complex_packed beta,
                          gsl_complex_packed_array C, int ldc)
{
  cblas_zsyr2k(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* TRMM */

void gsl_blas_raw_strmm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         float alpha,
                         const float A[], int lda,
                         float B[], int ldb)
{
  cblas_strmm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_dtrmm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         double alpha,
                         const double A[], int lda,
                         double B[], int ldb)
{
  cblas_dtrmm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ctrmm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float B, int ldb)
{
  cblas_ctrmm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ztrmm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array B, int ldb)
{
  cblas_ztrmm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


/* TRSM */

void gsl_blas_raw_strsm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         float alpha,
                         const float A[], int lda,
                         float B[], int ldb)
{
  cblas_strsm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_dtrsm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         double alpha,
                         const double A[], int lda,
                         double B[], int ldb)
{
  cblas_dtrsm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ctrsm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float B, int ldb)
{
  cblas_ctrsm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ztrsm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                         CBLAS_DIAG_t Diag,
                         size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array B, int ldb)
{
  cblas_ztrsm(CblasRowMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


/* HEMM */

void gsl_blas_raw_chemm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo,
                         size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_array_float B, int ldb,
                         const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
  cblas_chemm(CblasRowMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zhemm (CBLAS_SIDE_t Side,
                         CBLAS_UPLO_t Uplo,
                         size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array B, int ldb,
                         const gsl_const_complex_packed beta,
                         gsl_complex_packed_array C, int ldc)
{
  cblas_zhemm(CblasRowMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* HERK */

void gsl_blas_raw_cherk (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t Trans,
                         size_t N, size_t K,
                         float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
                         float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
  cblas_cherk(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_zherk (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t Trans,
                         size_t N, size_t K,
                         double alpha,
                         const gsl_const_complex_packed_array A, int lda,
                         double beta,
                         gsl_complex_packed_array C, int ldc)
{
  cblas_zherk(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


/* HER2K */

void gsl_blas_raw_cher2k (CBLAS_UPLO_t Uplo,
                          CBLAS_TRANSPOSE_t Trans,
                          size_t N, size_t K,
                          const gsl_const_complex_packed_float alpha,
                          const gsl_const_complex_packed_array_float A, int lda,
                          const gsl_const_complex_packed_array_float B, int ldb,
                          float beta,
                          gsl_complex_packed_array_float C, int ldc)
{
  cblas_cher2k(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zher2k (CBLAS_UPLO_t Uplo,
                          CBLAS_TRANSPOSE_t Trans,
                          size_t N, size_t K,
                          const gsl_const_complex_packed alpha,
                          const gsl_const_complex_packed_array A, int lda,
                          const gsl_const_complex_packed_array B, int ldb,
                          double beta,
                          gsl_complex_packed_array C, int ldc)
{
  cblas_zher2k(CblasRowMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
