/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_BLAS_H_
#define GSL_BLAS_H_

#include <gsl_matrix.h>
#include <gsl_vector.h>

/* !! IN PROGRESS: STILL NEEDS ALOT OF ARGUMENT MASSAGING 
   
   Must decide what to do with incX, incY args and with ld? args
   
   Also: Do we have a float complex, or just double?
 */



/* ========================================================================
 * Level 1
 * ========================================================================
 */

float gsl_blas_sdsdot (float alpha,
                       const gsl_vector_float * X, int incX,
                       const gsl_vector_float * Y, int incY);

double gsl_blas_dsdot (const gsl_vector_float X, int incX,
                       const gsl_vector_float Y, int incY);

float gsl_blas_sdot (const gsl_vector_float X, int incX,
                     const gsl_vector_float Y, int incY);

double gsl_blas_ddot (const gsl_vector X, int incX,
                      const gsl_vector Y, int incY);

/*
 * Functions having prefixes Z and C only
 */
void gsl_blas_cdotu (const gsl_vector_complex * X, int incX,
                     const gsl_vector_complex * Y, int incY,
                     gsl_vector_complex * dotu);

void gsl_blas_cdotc (const gsl_vector_complex * X, int incX,
                     const gsl_vector_complex * Y, int incY,
                     gsl_vector_complex * dotc);

void gsl_blas_zdotu (const gsl_vector_complex * X, int incX,
                     const gsl_vector_complex * Y, int incY,
                     gsl_vector_complex * dotu);

void gsl_blas_zdotc (const gsl_vector_complex * X, int incX,
                     const gsl_vector_complex * Y, int incY,
                     gsl_vector_complex * dotc);

/*
 * Functions having prefixes S D SC DZ
 */
float  gsl_blas_snrm2  (const gsl_vector_float * X, int incX);
float  gsl_blas_sasum  (const gsl_vector_float * X, int incX);
double gsl_blas_dnrm2  (const gsl_vector * X, int incX);
double gsl_blas_dasum  (const gsl_vector * X, int incX);
float  gsl_blas_scnrm2 (const gsl_vector_complex * X, int incX);
float  gsl_blas_scasum (const gsl_vector_complex * X, int incX);
double gsl_blas_dznrm2 (const gsl_vector_complex * X, int incX);
double gsl_blas_dzasum (const gsl_vector_complex * X, int incX);

/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX gsl_blas_isamax (const gsl_vector_float * X, int incX);
CBLAS_INDEX gsl_blas_idamax (const gsl_vector * X, int incX);
CBLAS_INDEX gsl_blas_icamax (const gsl_vector_complex * X, int incX);
CBLAS_INDEX gsl_blas_izamax (const gsl_vector_complex * X, int incX);


/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void gsl_blas_sswap (gsl_vector_float * X, int incX,
                     gsl_vector_float * Y, int incY);

void gsl_blas_scopy (const gsl_vector_float * X, int incX,
                     gsl_vector_float * Y, int incY);

void gsl_blas_saxpy (float alpha,
                     const gsl_vector_float * X, int incX,
                     gsl_vector_float * Y, int incY);

void gsl_blas_dswap (gsl_vector * X, int incX,
                     gsl_vector * Y, int incY);

void gsl_blas_dcopy (const gsl_vector * X, int incX,
                     gsl_vector * Y, int incY);

void gsl_blas_daxpy (double alpha,
                     const gsl_vector * X, int incX,
                     gsl_vector * Y, int incY);

void gsl_blas_cswap (gsl_vector_complex * X, int incX,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_ccopy (const gsl_vector_complex * X, int incX,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_caxpy (const gsl_vector_complex * alpha,
                     const gsl_vector_complex * X, int incX,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_zswap (gsl_vector_complex * X, int incX,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_zcopy (const gsl_vector_complex * X, int incX,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_zaxpy (const gsl_complex * alpha,
                     const gsl_vector_complex * X, int incX,
                     gsl_vector_complex * Y, int incY);

/*
 * Routines with S and D prefix only
 */
void gsl_blas_srotg (float a[], float b[], float c[], float s[]);

void gsl_blas_srotmg (float d1[], float d2[], float b1[], float b2, float P[]);

void gsl_blas_srot (gsl_vector_float * X, int incX,
                    gsl_vector_float * Y, int incY,
                    float c, float s);

void gsl_blas_srotm (gsl_vector_float * X, int incX,
                     gsl_vector_float * Y, int incY,
                     const float P[]);

void gsl_blas_drotg (double a[], double b[], double c[], double s[]);

void gsl_blas_drotmg (double d1[], double d2[], double b1[],
                      double b2, double P[]);

void gsl_blas_drot (gsl_vector * X, int incX,
                    gsl_vector * Y, int incY,
                    const double c, const double s);

void gsl_blas_drotm (gsl_vector * X, int incX,
                     gsl_vector * Y, int incY,
                     const double P[]);

/*
 * Routines with S D C Z CS and ZD prefixes
 */
void gsl_blas_sscal  (float  alpha, gsl_vector_float  * X, int incX);
void gsl_blas_dscal  (double alpha, gsl_vector * X, int incX);
void gsl_blas_cscal  (const gsl_complex * alpha, gsl_vector_complex * X, int incX);
void gsl_blas_zscal  (const gsl_complex * alpha, gsl_vector_complex * X, int incX);
void gsl_blas_csscal (float  alpha, gsl_vector_complex * X, int incX);
void gsl_blas_zdscal (double alpha, gsl_vector_complex * X, int incX);



/* ===========================================================================
 * Level 2
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void gsl_blas_sgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int M, int N,
                     float alpha,
                     const gsl_matrix_float * A, int lda,
                     const gsl_vector_float * X, int incX,
                     float beta,
                     gsl_vector_float * Y, int incY);

void gsl_blas_sgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int M, int N, int KL, int KU,
		     float alpha,
                     const gsl_matrix_float * A, int lda,
                     const gsl_vector_float * X, int incX,
                     float beta,
                     gsl_vector_float * Y, int incY);

void gsl_blas_strmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const gsl_matrix_float * A, int lda,
                     gsl_vector_float * X, int incX);

void gsl_blas_stbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const gsl_matrix_float * A, int lda,
                     gsl_vector_float * X, int incX);

void gsl_blas_stpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const float Ap[],
                     gsl_vector_float * X, int incX);

void gsl_blas_strsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const float A[], int lda,
                     gsl_vector_float * X, int incX);

void gsl_blas_stbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const float A[], int lda,
                     gsl_vector_float * X, int incX);

void gsl_blas_stpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const float Ap[],
                     gsl_vector_float * X, int incX);

void gsl_blas_dgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int M, int N,
                     double alpha,
                     const double A[], int lda,
                     const gsl_vector * X, int incX,
                     double beta,
                     gsl_vector * Y, int incY);

void gsl_blas_dgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int M, int N, int KL, int KU,
                     double alpha,
                     const double A[], int lda,
                     const gsl_vector * X, int incX,
                     double beta,
                     gsl_vector * Y, int incY);

void gsl_blas_dtrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const double A[], int lda,
                     gsl_vector * X, int incX);

void gsl_blas_dtbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const double A[], int lda,
                     gsl_vector * X, int incX);

void gsl_blas_dtpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const double Ap[],
                     gsl_vector * X, int incX);

void gsl_blas_dtrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const double A[], int lda,
                     gsl_vector * X, int incX);

void gsl_blas_dtbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const double A[], int lda,
                     gsl_vector * X, int incX);

void gsl_blas_dtpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const double Ap[],
                     gsl_vector * X, int incX);

void gsl_blas_cgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void * A, int lda,
                     const gsl_vector_complex * X, int incX,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_cgbmv (CBLAS_ORDER order,
                     CBLAS_TRANSPOSE TransA,
                     int M, int N, int KL, int KU,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const gsl_vector_complex * X, int incX,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_ctrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *A, int lda,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ctbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const void *A, int lda,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ctpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *Ap,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ctrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *A, int lda,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ctbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const void *A, int lda,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ctpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *Ap,
                     gsl_vector_complex * X, int incX);

void gsl_blas_zgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const gsl_vector_complex * X, int incX,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_zgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int M, int N, int KL, int KU,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const gsl_vector_complex * X, int incX,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_ztrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *A, int lda,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ztbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const void *A, int lda,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ztpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *Ap,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ztrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *A, int lda,
                     gsl_vector_complex *X, int incX);

void gsl_blas_ztbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N, int K,
                     const void *A, int lda,
                     gsl_vector_complex * X, int incX);

void gsl_blas_ztpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int N,
                     const void *Ap,
                     gsl_vector_complex * X, int incX);

/*
 * Routines with S and D prefixes only
 */
void gsl_blas_ssymv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     float alpha,
                     const float A[], int lda,
                     const gsl_vector_float * X, int incX,
                     float beta,
                     gsl_vector_float * Y, int incY);

void gsl_blas_ssbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N, int K,
                     float alpha,
                     const float A[], int lda,
                     const gsl_vector_float * X, int incX,
                     float beta,
                     gsl_vector_float * Y, int incY);

void gsl_blas_sspmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     float alpha,
                     const float Ap[],
                     const gsl_vector_float * X, int incX,
                     float beta,
                     gsl_vector_float * Y, int incY);

void gsl_blas_sger (CBLAS_ORDER order,
                    int M, int N,
                    float alpha,
                    const gsl_vector_float * X, int incX,
                    const gsl_vector_float * Y, int incY,
                    float A[], int lda);

void gsl_blas_ssyr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    float alpha,
                    const gsl_vector_float * X, int incX,
                    float A[], int lda);

void gsl_blas_sspr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    float alpha,
                    const gsl_vector_float * X, int incX,
                    float Ap[]);

void gsl_blas_ssyr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     float alpha,
                     const gsl_vector_float * X, int incX,
                     const gsl_vector_float * Y, int incY,
                     float A[], int lda);

void gsl_blas_sspr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     float alpha,
                     const gsl_vector_float * X, int incX,
                     const gsl_vector_float * Y, int incY,
                     gsl_vector_float * A);

void gsl_blas_dsymv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     double alpha,
                     const double A[], int lda,
                     const gsl_vector * X, int incX,
                     double beta,
                     gsl_vector * Y, int incY);

void gsl_blas_dsbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N, int K,
                     double alpha,
                     const double A[], int lda,
                     const gsl_vector * X, int incX,
                     double beta,
                     gsl_vector * Y, int incY);

void gsl_blas_dspmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     double alpha,
                     const double Ap[],
                     const gsl_vector * X, int incX,
                     double beta,
                     gsl_vector * Y, int incY);

void gsl_blas_dger (CBLAS_ORDER order,
                    int M, int N,
                    double alpha,
                    const gsl_vector * X, int incX,
                    const gsl_vector * Y, int incY,
                    double A[], int lda);

void gsl_blas_dsyr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    double alpha,
                    const gsl_vector * X, int incX,
                    double A[], int lda);

void gsl_blas_dspr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    double alpha,
                    const gsl_vector * X, int incX,
                    double Ap[]);

void gsl_blas_dsyr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     double alpha,
                     const gsl_vector * X, int incX,
                     const gsl_vector * Y, int incY,
                     double A[], int lda);

void gsl_blas_dspr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     double alpha,
                     const gsl_vector * X, int incX,
                     const gsl_vector * Y, int incY,
                     double A[]);

/*
 * Routines with C and Z prefixes only
 */
void gsl_blas_chemv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const gsl_vector_complex * X, int incX,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_chbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N, int K,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const gsl_vector_complex * X, int incX,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y, int incY);

void gsl_blas_chpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
		     const void *Ap,
                     const void *X, int incX,
                     const gsl_complex * beta,
                     void *Y, int incY);

void gsl_blas_cgeru (CBLAS_ORDER order,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *A, int lda);

void gsl_blas_cgerc (CBLAS_ORDER order,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *A, int lda);

void gsl_blas_cher (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    float alpha,
                    const void *X, int incX,
                    void *A, int lda);

void gsl_blas_chpr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    float alpha,
                    const void *X, int incX,
                    void *A);

void gsl_blas_cher2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *A, int lda);

void gsl_blas_chpr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *Ap);

void gsl_blas_zhemv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *X, int incX,
                     const gsl_complex * beta,
                     void *Y, int incY);

void gsl_blas_zhbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N, int K,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *X, int incX,
                     const gsl_complex * beta,
                     void *Y, int incY);

void gsl_blas_zhpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
                     const void *Ap,
                     const void *X, int incX,
                     const gsl_complex * beta,
                     void *Y, int incY);

void gsl_blas_zgeru (CBLAS_ORDER order,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *A, int lda);

void gsl_blas_zgerc (CBLAS_ORDER order,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *A, int lda);

void gsl_blas_zher (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    double alpha,
                    const void *X, int incX,
                    void *A, int lda);

void gsl_blas_zhpr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    int N,
                    double alpha,
                    const void *X, int incX,
                    void *A);

void gsl_blas_zher2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *A, int lda);

void gsl_blas_zhpr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int N,
                     const gsl_complex * alpha,
                     const void *X, int incX,
                     const void *Y, int incY,
                     void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void gsl_blas_sgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int M, int N, int K,
                     float alpha,
                     const float A[], int lda,
                     const float B[], int ldb,
                     float beta,
                     float C[], int ldc);

void gsl_blas_ssymm (CBLAS_ORDER Order, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                     int M, int N,
                     float alpha,
                     const float A[], int lda,
                     const float B[], int ldb,
                     float beta,
                     float C[], int ldc);

void gsl_blas_ssyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                     int N, int K,
                     float alpha,
                     const float A[], int lda,
                     float beta,
                     float C[], int ldc);

void gsl_blas_ssyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                      int N, int K,
                      float alpha,
                      const float A[], int lda,
                      const float B[], int ldb,
                      float beta,
                      float C[], int ldc);

void gsl_blas_strmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     float alpha,
                     const float A[], int lda,
                     float B[], int ldb);

void gsl_blas_strsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     float alpha,
                     const float A[], int lda,
                     float B[], int ldb);

void gsl_blas_dgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int M, int N, int K,
                     double alpha,
                     const double A[], int lda,
                     const double B[], int ldb,
                     double beta,
                     double C[], int ldc);

void gsl_blas_dsymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     int M, int N,
                     double alpha,
                     const double A[], int lda,
                     const double B[], int ldb,
                     double beta,
                     double C[], int ldc);

void gsl_blas_dsyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int N, int K,
                     double alpha,
                     const double A[], int lda,
                     double beta,
                     double C[], int ldc);

void gsl_blas_dsyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int N, int K,
                      double alpha,
                      const double A[], int lda,
                      const double B[], int ldb,
                      double beta,
                      double C[], int ldc);

void gsl_blas_dtrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     double alpha,
                     const double A[], int lda,
                     double B[], int ldb);

void gsl_blas_dtrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     double alpha,
                     const double A[], int lda,
                     double B[], int ldb);

void gsl_blas_cgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int M, int N, int K,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *B, int ldb,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_csymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *B, int ldb,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_csyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int N, int K,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_csyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int N, int K,
                      const gsl_complex * alpha,
                      const void *A, int lda,
                      const void *B, int ldb,
                      const gsl_complex * beta,
                      void *C, int ldc);

void gsl_blas_ctrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     void *B, int ldb);

void gsl_blas_ctrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     void *B, int ldb);

void gsl_blas_zgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int M, int N, int K,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *B, int ldb,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_zsymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *B, int ldb,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_zsyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int N, int K,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_zsyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int N, int K,
                      const gsl_complex * alpha,
                      const void *A, int lda,
                      const void *B, int ldb,
                      const gsl_complex * beta,
                      void *C, int ldc);

void gsl_blas_ztrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     void *B, int ldb);

void gsl_blas_ztrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     void *B, int ldb);

/*
 * Routines with prefixes C and Z only
 */
void gsl_blas_chemm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *B, int ldb,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_cherk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int N, int K,
                     float alpha,
                     const void *A, int lda,
                     float beta,
                     void *C, int ldc);

void gsl_blas_cher2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int N, int K,
                      const gsl_complex * alpha,
                      const void *A, int lda,
                      const void *B, int ldb,
                      float beta,
                      void *C, int ldc);

void gsl_blas_zhemm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     int M, int N,
                     const gsl_complex * alpha,
                     const void *A, int lda,
                     const void *B, int ldb,
                     const gsl_complex * beta,
                     void *C, int ldc);

void gsl_blas_zherk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int N, int K,
                     double alpha,
                     const void *A, int lda,
                     double beta,
                     void *C, int ldc);

void gsl_blas_zher2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int N, int K,
                      const gsl_complex * alpha,
                      const void *A, int lda,
                      const void *B, int ldb,
                      double beta,
                      void *C, int ldc);



#endif /* !GSL_BLAS_H_ */
