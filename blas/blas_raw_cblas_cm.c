/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Implementation of gsl_blas_raw_cm interface which
 * defers to a CBLAS conformant interface. This includes
 * only the level 2 and level 3 functions, since the level 1
 * functions are insensitive to the storage scheme.
 */
#include "gsl_blas_raw_cm.h"


/*
 * ===========================================================================
 * level 2 BLAS
 * ===========================================================================
 */

 
/* GEMV */

void gsl_blas_raw_sgemv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY)
{
  cblas_sgemv(CblasColMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dgemv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY)
{
  cblas_dgemv(CblasColMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_cgemv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_cgemv(CblasColMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zgemv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_zgemv(CblasColMajor, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
}


/* GBMV */

void gsl_blas_raw_sgbmv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY)
{
  cblas_sgbmv(CblasColMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dgbmv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY)
{
  cblas_dgbmv(CblasColMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_cgbmv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_cgbmv(CblasColMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zgbmv_cm (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_zgbmv(CblasColMajor, TransA, M, N, KL, KU, alpha, A, lda, X, incX, beta, Y, incY);
}


/* TRMV */

void gsl_blas_raw_strmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float A[], int lda,
                            float X[], int incX)
{
  cblas_strmv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_dtrmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double A[], int lda,
                            double X[], int incX)
{
  cblas_dtrmv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ctrmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ctrmv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ztrmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ztrmv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


/* TBMV */

void gsl_blas_raw_stbmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const float A[], int lda,
                            float X[], int incX)
{
  cblas_stbmv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_dtbmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const double A[], int lda,
                            double X[], int incX)
{
  cblas_dtbmv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ctbmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ctbmv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ztbmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ztbmv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


/* TPMV */

void gsl_blas_raw_stpmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float Ap[],
                            float X[], int incX)
{
  cblas_stpmv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_dtpmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double Ap[],
                            double X[], int incX)
{
  cblas_dtpmv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ctpmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX)
{
  cblas_ctpmv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ztpmv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX)
{
  cblas_ztpmv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


/* TRSV */

void gsl_blas_raw_strsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float A[], int lda,
                            float X[], int incX)
{
  cblas_strsv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_dtrsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double A[], int lda,
                            double X[], int incX)
{
  cblas_dtrsv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ctrsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ctrsv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


void gsl_blas_raw_ztrsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ztrsv(CblasColMajor, Uplo, TransA, Diag, N, A, lda, X, incX);
}


/* TBSV */

void gsl_blas_raw_stbsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const float A[], int lda,
                            float X[], int incX)
{
  cblas_stbsv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_dtbsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const double A[], int lda,
                            double X[], int incX)
{
  cblas_dtbsv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ctbsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ctbsv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


void gsl_blas_raw_ztbsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX)
{
  cblas_ztbsv(CblasColMajor, Uplo, TransA, Diag, N, K, A, lda, X, incX);
}


/* TPSV */

void gsl_blas_raw_stpsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float Ap[],
                            float X[], int incX)
{
  cblas_stpsv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_dtpsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double Ap[],
                            double X[], int incX)
{
  cblas_dtpsv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ctpsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX)
{
  cblas_ctpsv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


void gsl_blas_raw_ztpsv_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX)
{
  cblas_ztpsv(CblasColMajor, Uplo, TransA, Diag, N, Ap, X, incX);
}


/* SYMV */

void gsl_blas_raw_ssymv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY)
{
  cblas_ssymv(CblasColMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dsymv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY)
{
  cblas_dsymv(CblasColMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


/* SBMV */

void gsl_blas_raw_ssbmv_cm (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            float alpha,
                            const float A[], int lda, 
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY)
{
  cblas_ssbmv(CblasColMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dsbmv_cm (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY)
{
  cblas_dsbmv(CblasColMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


/* SPMV */

void gsl_blas_raw_sspmv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                            const float Ap[],
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY)
{
  cblas_sspmv(CblasColMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


void gsl_blas_raw_dspmv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double Ap[],
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY)
{
  cblas_dspmv(CblasColMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


/* GER */

void gsl_blas_raw_sger_cm (size_t M, size_t N,
                           float alpha,
                           const float X[], int incX,
                           const float Y[], int incY,
                           float A[], int lda)
{
  cblas_sger(CblasColMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_dger_cm (size_t M, size_t N,
                           double alpha,
                           const double X[], int incX,
                           const double Y[], int incY,
                           double A[], int lda)
{
  cblas_dger(CblasColMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


/* SYR */

void gsl_blas_raw_ssyr_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const float X[], int incX,
                           float A[], int lda)
{
  cblas_ssyr(CblasColMajor, Uplo, N,alpha, X, incX, A, lda);
}


void gsl_blas_raw_dsyr_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const double X[], int incX,
                           double A[], int lda)
{
  cblas_dsyr(CblasColMajor, Uplo, N, alpha, X, incX, A, lda);
}


/* SPR */

void gsl_blas_raw_sspr_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const float X[], int incX,
                           float Ap[])
{
  cblas_sspr(CblasColMajor, Uplo, N, alpha, X, incX, Ap);
}


void gsl_blas_raw_dspr_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const double X[], int incX,
                           double Ap[])
{
  cblas_dspr(CblasColMajor, Uplo, N, alpha, X, incX, Ap);
}


/* SYR2 */

void gsl_blas_raw_ssyr2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                            const float X[], int incX, 
                            const float Y[], int incY,
                            float A[], int lda)
{
  cblas_ssyr2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_dsyr2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double X[], int incX,
                            const double Y[], int incY,
                            double A[], int lda)
{
  cblas_dsyr2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


/* SPR2 */

void gsl_blas_raw_sspr2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                            const float X[], int incX,
                            const float Y[], int incY,
                            float A[])
{
  cblas_sspr2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, A);
}


void gsl_blas_raw_dspr2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double X[], int incX,
                            const double Y[], int incY,
                            double A[])
{
  cblas_sspr2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, A);
}


/* HEMV */

void gsl_blas_raw_chemv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_chemv(CblasColMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zhemv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_zhemv(CblasColMajor, Uplo, N, alpha, A, lda, X, incX, beta, Y, incY);
}


/* HBMV */

void gsl_blas_raw_chbmv_cm (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_chbmv(CblasColMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zhbmv_cm (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_zhbmv(CblasColMajor, Uplo, N, K, alpha, A, lda, X, incX, beta, Y, incY);
}


/* HPMV */

void gsl_blas_raw_chpmv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha, const void * Ap,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_chpmv(CblasColMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


void gsl_blas_raw_zhpmv_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * Ap,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY)
{
  cblas_zhpmv(CblasColMajor, Uplo, N, alpha, Ap, X, incX, beta, Y, incY);
}


/* GERU */

void gsl_blas_raw_cgeru_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda)
{
  cblas_cgeru(CblasColMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_zgeru_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda)
{
  cblas_zgeru(CblasColMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


/* GERC */

void gsl_blas_raw_cgerc_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda)
{
  cblas_cgerc(CblasColMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_zgerc_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda)
{
  cblas_zgerc(CblasColMajor, M, N, alpha, X, incX, Y, incY, A, lda);
}


/* HER */

void gsl_blas_raw_cher_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const void * X, int incX,
                           void * A, int lda)
{
  cblas_cher(CblasColMajor, Uplo, N, X, incX, A, lda);
}


void gsl_blas_raw_zher_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const void * X, int incX,
                           void * A, int lda)
{
  cblas_zher(CblasColMajor, Uplo, N, X, incX, A, lda);
}


/* HPR */

void gsl_blas_raw_chpr_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const void * X, int incX,
                           void * A)
{
  cblas_chpr(CblasColMajor, Uplo, N, alpha, X, incX, A);
}


void gsl_blas_raw_zhpr_cm (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const void * X, int incX,
                           void * A)
{
  cblas_zhpr(CblasColMajor, Uplo, N, alpha, X, incX, A);
}


/* HER2 */

void gsl_blas_raw_cher2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda)
{
  cblas_cher2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


void gsl_blas_raw_zher2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda)
{
  cblas_zher2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, A, lda);
}


/* HPR2 */

void gsl_blas_raw_chpr2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * Ap)
{
  cblas_chpr2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, Ap);
}


void gsl_blas_raw_zhpr2_cm (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * Ap)
{
  cblas_zhpr2(CblasColMajor, Uplo, N, alpha, X, incX, Y, incY, Ap);
}



/*
 * ===========================================================================
 * level 3 BLAS
 * ===========================================================================
 */


/* GEMM */

void gsl_blas_raw_sgemm_cm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            float alpha,
                            const float A[], int lda,
                            const float B[], int ldb,
                            float beta,
                            float C[], int ldc)
{
  cblas_sgemm(CblasColMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_dgemm_cm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            const double B[], int ldb,
                            double beta,
                            double C[], int ldc)
{
  cblas_dgemm(CblasColMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_cgemm_cm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda, 
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_cgemm(CblasColMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zgemm_cm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_zgemm(CblasColMajor, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* SYMM */

void gsl_blas_raw_ssymm_cm (CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float B[], int ldb,
                            float beta,
                            float C[], int ldc)
{
  cblas_ssymm(CblasColMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_dsymm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double B[], int ldb,
                            double beta,
                            double C[], int ldc)
{
  cblas_dsymm(CblasColMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);

}


void gsl_blas_raw_csymm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_csymm(CblasColMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zsymm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_zsymm(CblasColMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* SYRK */

void gsl_blas_raw_ssyrk_cm (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            float alpha,
                            const float A[], int lda,
                            float beta,
                            float C[], int ldc)
{
  cblas_ssyrk(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_dsyrk_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            double beta,
                            double C[], int ldc)
{
  cblas_dsyrk(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_csyrk_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_csyrk(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_zsyrk_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_zsyrk(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


/* SYR2K */

void gsl_blas_raw_ssyr2k_cm (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             float alpha,
                             const float A[], int lda,
                             const float B[], int ldb,
                             float beta,
                             float C[], int ldc)
{
  cblas_ssyr2k(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_dsyr2k_cm (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             double alpha,
                             const double A[], int lda,
                             const double B[], int ldb,
                             double beta,
                             double C[], int ldc)
{
  cblas_dsyr2k(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_csyr2k_cm (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             const void * beta,
                             void * C, int ldc)
{
  cblas_csyr2k(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zsyr2k_cm (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             const void * beta,
                             void * C, int ldc)
{
  cblas_zsyr2k(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* TRMM */

void gsl_blas_raw_strmm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            float B[], int ldb)
{
  cblas_strmm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_dtrmm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            double B[], int ldb)
{
  cblas_dtrmm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ctrmm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb)
{
  cblas_ctrmm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ztrmm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb)
{
  cblas_ztrmm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


/* TRSM */

void gsl_blas_raw_strsm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            float B[], int ldb)
{
  cblas_strsm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_dtrsm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            double B[], int ldb)
{
  cblas_dtrsm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ctrsm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb)
{
  cblas_ctrsm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


void gsl_blas_raw_ztrsm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb)
{
  cblas_ztrsm(CblasColMajor, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
}


/* HEMM */

void gsl_blas_raw_chemm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_chemm(CblasColMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zhemm_cm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc)
{
  cblas_zhemm(CblasColMajor, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* HERK */

void gsl_blas_raw_cherk_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            float alpha,
                            const void * A, int lda,
                            float beta,
                            void * C, int ldc)
{
  cblas_cherk(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


void gsl_blas_raw_zherk_cm (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            double alpha,
                            const void * A, int lda,
                            double beta,
                            void * C, int ldc)
{
  cblas_zherk(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
}


/* HER2K */

void gsl_blas_raw_cher2k_cm (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             float beta,
                             void * C, int ldc)
{
  cblas_cher2k(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}


void gsl_blas_raw_zher2k_cm (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             double beta,
                             void * C, int ldc)
{
  cblas_zher2k(CblasColMajor, Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
