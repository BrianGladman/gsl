/*
 * Author:  G. Jungman  and  B. Gough
 * RCS:     $Id$
 */
#ifndef GSL_BLAS_NATIVE_H_
#define GSL_BLAS_NATIVE_H_

#include "gsl_blas_types.h"


/*
 * ===========================================================================
 * level 1 BLAS functions
 * ===========================================================================
 */

float  gsl_blas_native_sdsdot (size_t N,
                               float alpha,
                               const float X[], int incX,
                               const float Y[], int incY);

double gsl_blas_native_dsdot (size_t N,
                              const float X[], int incX,
                              const float Y[], int incY);

float  gsl_blas_native_sdot (size_t N,
                             const float X[], int incX,
                             const float Y[], int incY);
                          
double gsl_blas_native_ddot (size_t N,
                             const double X[], int incX,
                             const double Y[], int incY);

void gsl_blas_native_cdotu (size_t N,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * dotu);

void gsl_blas_native_cdotc (size_t N,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * dotc);

void gsl_blas_native_zdotu (size_t N,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * dotu);

void gsl_blas_native_zdotc (size_t N,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * dotc);

float  gsl_blas_native_snrm2  (size_t N, const float  X[], int incX);
double gsl_blas_native_dnrm2  (size_t N, const double X[], int incX);
float  gsl_blas_native_scnrm2 (size_t N, const void * X, int incX);
double gsl_blas_native_dznrm2 (size_t N, const void * X, int incX);
float  gsl_blas_native_sasum  (size_t N, const float  X[], int incX);
double gsl_blas_native_dasum  (size_t N, const double X[], int incX);
float  gsl_blas_native_scasum (size_t N, const void * X, int incX);
double gsl_blas_native_dzasum (size_t N, const void * X, int incX);

CBLAS_INDEX gsl_blas_native_isamax (size_t N, const float  X[], int incX);
CBLAS_INDEX gsl_blas_native_idamax (size_t N, const double X[], int incX);
CBLAS_INDEX gsl_blas_native_icamax (size_t N, const void * X, int incX);
CBLAS_INDEX gsl_blas_native_izamax (size_t N, const void * X, int incX);


/*
 * ===========================================================================
 * level 1 BLAS routines
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void gsl_blas_native_sswap (size_t N,
                            float X[], int incX,
                            float Y[], int incY);

void gsl_blas_native_dswap (size_t N,
                            double X[], int incX,
                            double Y[], int incY);

void gsl_blas_native_cswap (size_t N,
                            void * X, int incX,
                            void * Y, int incY);

void gsl_blas_native_zswap (size_t N,
                            void * X, int incX,
                            void * Y, int incY);

void gsl_blas_native_scopy (size_t N,
                            const float X[], int incX,
                            float Y[], int incY);

void gsl_blas_native_dcopy (size_t N,
                            const double X[], int incX,
                            double Y[], int incY);

void gsl_blas_native_ccopy (size_t N,
                            const void * X, int incX,
                            void * Y, int incY);

void gsl_blas_native_zcopy (size_t N,
                            const void * X, int incX,
                            void * Y, int incY);

void gsl_blas_native_saxpy (size_t N,
                            float alpha,
                            const float X[], int incX,
                            float Y[], int incY);


void gsl_blas_native_daxpy (size_t N,
                            double alpha,
                            const double X[], int incX, 
                            double Y[], int incY);


void gsl_blas_native_caxpy (size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            void * Y, int incY);

void gsl_blas_native_zaxpy (size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            void * Y, int incY);

/*
 * Routines with S and D prefix only
 */
void gsl_blas_native_srotg (float a[], float b[], float c[], float s[]);
void gsl_blas_native_drotg (double a[], double b[], double c[], double s[]);

void gsl_blas_native_srotmg (float d1[], float d2[], float b1[], float b2, float P[]);
void gsl_blas_native_drotmg (double d1[], double d2[], double b1[],
                             double b2, double P[]);


void gsl_blas_native_srot (size_t N,
                           float X[], int incX,
                           float Y[], int incY,
                           float c, float s);

void gsl_blas_native_drot (size_t N,
                           double X[], int incX,
                           double Y[], int incY,
                           const double c, const double s);

void gsl_blas_native_srotm (size_t N,
                            float X[], int incX,
                            float Y[], int incY,
                            const float P[]);

void gsl_blas_native_drotm (size_t N,
                            double X[], int incX,
                            double Y[], int incY,
                            const double P[]);

/*
 * Routines with S D C Z CS and ZD prefixes
 */
void gsl_blas_native_sscal  (size_t N, float  alpha, float  X[], int incX);
void gsl_blas_native_dscal  (size_t N, double alpha, double X[], int incX);
void gsl_blas_native_cscal  (size_t N, const void * alpha, void * X, int incX);
void gsl_blas_native_zscal  (size_t N, const void * alpha, void * X, int incX);
void gsl_blas_native_csscal (size_t N, float  alpha, void * X, int incX);
void gsl_blas_native_zdscal (size_t N, double alpha, void * X, int incX);


/*
 * ===========================================================================
 * level 2 BLAS
 * ===========================================================================
 */

 
/* GEMV */

void gsl_blas_native_sgemv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_native_dgemv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

void gsl_blas_native_cgemv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_native_zgemv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

/* GBMV */

void gsl_blas_native_sgbmv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_native_dgbmv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

void gsl_blas_native_cgbmv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_native_zgbmv (CBLAS_TRANSPOSE TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* TRMV */

void gsl_blas_native_strmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_native_dtrmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_native_ctrmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_native_ztrmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);


/* TBMV */

void gsl_blas_native_stbmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_native_dtbmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_native_ctbmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_native_ztbmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);


/* TPMV */

void gsl_blas_native_stpmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float Ap[],
                            float X[], int incX);

void gsl_blas_native_dtpmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double Ap[],
                            double X[], int incX);

void gsl_blas_native_ctpmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);

void gsl_blas_native_ztpmv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);

/* TRSV */

void gsl_blas_native_strsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_native_dtrsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_native_ctrsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_native_ztrsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);


/* TBSV */

void gsl_blas_native_stbsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_native_dtbsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_native_ctbsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_native_ztbsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);


/* TPSV */

void gsl_blas_native_stpsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const float Ap[],
                            float X[], int incX);

void gsl_blas_native_dtpsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const double Ap[],
                            double X[], int incX);

void gsl_blas_native_ctpsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);

void gsl_blas_native_ztpsv (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);


/* SYMV */

void gsl_blas_native_ssymv (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                               const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_native_dsymv (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);


/* SBMV */

void gsl_blas_native_ssbmv (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            float alpha,
                            const float A[], int lda, 
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_native_dsbmv (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

/* SPMV */

void gsl_blas_native_sspmv (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                            const float Ap[],
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_native_dspmv (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double Ap[],
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

/* GER */

void gsl_blas_native_sger (CBLAS_ORDER order,
                           size_t M, size_t N,
                           float alpha,
                           const float X[], int incX,
                           const float Y[], int incY,
                           float A[], int lda);

void gsl_blas_native_dger (CBLAS_ORDER order,
                           size_t M, size_t N,
                           double alpha,
                           const double X[], int incX,
                           const double Y[], int incY,
                           double A[], int lda);


/* SYR */

void gsl_blas_native_ssyr (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const float X[], int incX,
                           float A[], int lda);

void gsl_blas_native_dsyr (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const double X[], int incX,
                           double A[], int lda);


/* SPR */

void gsl_blas_native_sspr (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const float X[], int incX,
                           float Ap[]);

void gsl_blas_native_dspr (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const double X[], int incX,
                           double Ap[]);


/* SYR2 */

void gsl_blas_native_ssyr2 (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                            const float X[], int incX, 
                            const float Y[], int incY,
                            float A[], int lda);

void gsl_blas_native_dsyr2 (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double X[], int incX,
                            const double Y[], int incY,
                            double A[], int lda);


/* SPR2 */

void gsl_blas_native_sspr2 (CBLAS_UPLO Uplo,
                            size_t N,
                            float alpha,
                            const float X[], int incX,
                            const float Y[], int incY,
                            float A[]);

void gsl_blas_native_dspr2 (CBLAS_UPLO Uplo,
                            size_t N,
                            double alpha,
                            const double X[], int incX,
                            const double Y[], int incY,
                            double A[]);


/* HEMV */

void gsl_blas_native_chemv (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_native_zhemv (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* HBMV */

void gsl_blas_native_chbmv (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_native_zhbmv (CBLAS_UPLO Uplo,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* HPMV */

void gsl_blas_native_chpmv (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha, const void * Ap,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_native_zhpmv (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * Ap,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* GERU */

void gsl_blas_native_cgeru (CBLAS_ORDER order,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

void gsl_blas_native_zgeru (CBLAS_ORDER order,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);


/* GERC */

void gsl_blas_native_cgerc (CBLAS_ORDER order,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

void gsl_blas_native_zgerc (CBLAS_ORDER order,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

/* HER */

void gsl_blas_native_cher (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const void * X, int incX,
                           void * A, int lda);

void gsl_blas_native_zher (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const void * X, int incX,
                           void * A, int lda);


/* HPR */

void gsl_blas_native_chpr (CBLAS_UPLO Uplo,
                           size_t N,
                           float alpha,
                           const void * X, int incX,
                           void * A);

void gsl_blas_native_zhpr (CBLAS_UPLO Uplo,
                           size_t N,
                           double alpha,
                           const void * X, int incX,
                           void * A);


/* HER2 */

void gsl_blas_native_cher2 (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

void gsl_blas_native_zher2 (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);


/* HPR2 */

void gsl_blas_native_chpr2 (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * Ap);

void gsl_blas_native_zhpr2 (CBLAS_UPLO Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * Ap);


/*
 * ===========================================================================
 * level 3 BLAS
 * ===========================================================================
 */

/* GEMM */

void gsl_blas_native_sgemm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            float alpha,
                            const float A[], int lda,
                            const float B[], int ldb,
                            float beta,
                            float C[], int ldc);

void gsl_blas_native_dgemm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            const double B[], int ldb,
                            double beta,
                            double C[], int ldc);

void gsl_blas_native_cgemm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda, 
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_native_zgemm (CBLAS_TRANSPOSE TransA,
                            CBLAS_TRANSPOSE TransB,
                            size_t M, size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);


/* SYMM */

void gsl_blas_native_ssymm (CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float B[], int ldb,
                            float beta,
                            float C[], int ldc);

void gsl_blas_native_dsymm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double B[], int ldb,
                            double beta,
                            double C[], int ldc);

void gsl_blas_native_csymm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_native_zsymm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);


/* SYRK */

void gsl_blas_native_ssyrk (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            float alpha,
                            const float A[], int lda,
                            float beta,
                            float C[], int ldc);

void gsl_blas_native_dsyrk (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            double beta,
                            double C[], int ldc);

void gsl_blas_native_csyrk (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_native_zsyrk (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * beta,
                            void * C, int ldc);


/* SYR2K */

void gsl_blas_native_ssyr2k (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             float alpha,
                             const float A[], int lda,
                             const float B[], int ldb,
                             float beta,
                             float C[], int ldc);

void gsl_blas_native_dsyr2k (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             double alpha,
                             const double A[], int lda,
                             const double B[], int ldb,
                             double beta,
                             double C[], int ldc);

void gsl_blas_native_csyr2k (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             const void * beta,
                             void * C, int ldc);

void gsl_blas_native_zsyr2k (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             const void * beta,
                             void * C, int ldc);


/* TRMM */

void gsl_blas_native_strmm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            float B[], int ldb);

void gsl_blas_native_dtrmm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            double B[], int ldb);

void gsl_blas_native_ctrmm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);

void gsl_blas_native_ztrmm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);


/* TRSM */

void gsl_blas_native_strsm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            float B[], int ldb);

void gsl_blas_native_dtrsm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            double B[], int ldb);

void gsl_blas_native_ctrsm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);

void gsl_blas_native_ztrsm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                            CBLAS_DIAG Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);


/* HEMM */

void gsl_blas_native_chemm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_native_zhemm (CBLAS_SIDE Side,
                            CBLAS_UPLO Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);


/* HERK */

void gsl_blas_native_cherk (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            float alpha,
                            const void * A, int lda,
                            float beta,
                            void * C, int ldc);

void gsl_blas_native_zherk (CBLAS_UPLO Uplo,
                            CBLAS_TRANSPOSE Trans,
                            size_t N, size_t K,
                            double alpha,
                            const void * A, int lda,
                            double beta,
                            void * C, int ldc);


/* HER2K */

void gsl_blas_native_cher2k (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             float beta,
                             void * C, int ldc);


void gsl_blas_native_zher2k (CBLAS_UPLO Uplo,
                             CBLAS_TRANSPOSE Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             double beta,
                             void * C, int ldc);


#endif  /* !GSL_BLAS_NATIVE_H_ */
