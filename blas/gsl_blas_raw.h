/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef GSL_BLAS_RAW_H_
#define GSL_BLAS_RAW_H_

#include <gsl_blas_types.h>


/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */

float  gsl_blas_raw_sdsdot (int N,
                            float alpha,
                            const float X[], int incX,
                            const float Y[], int incY);

double gsl_blas_raw_dsdot (int N,
                           const float X[], int incX,
                           const float Y[], int incY);

float  gsl_blas_raw_sdot (int N,
                          const float X[], int incX,
                          const float Y[], int incY);
                          
double gsl_blas_raw_ddot (int N,
                          const double X[], int incX,
                          const double Y[], int incY);

/*
 * Functions having prefixes Z and C only
 */
void gsl_blas_raw_cdotu (int N,
                         const void * X, int incX,
                         const void * Y, int incY,
		         void * dotu);

void gsl_blas_raw_cdotc (int N,
                         const void * X, int incX,
                         const void * Y, int incY,
                         void * dotc);

void gsl_blas_raw_zdotu (int N,
                         const void * X, int incX,
                         const void * Y, int incY,
                         void * dotu);

void gsl_blas_raw_zdotc (int N,
                         const void * X, int incX,
                         const void * Y, int incY,
                         void * dotc);

/*
 * Functions having prefixes S D SC DZ
 */
float  gsl_blas_raw_snrm2  (int N, const float X[], int incX);
float  gsl_blas_raw_sasum  (int N, const float X[], int incX);
double gsl_blas_raw_dnrm2  (int N, const double X[], int incX);
double gsl_blas_raw_dasum  (int N, const double X[], int incX);
float  gsl_blas_raw_scnrm2 (int N, const void * X, int incX);
float  gsl_blas_raw_scasum (int N, const void * X, int incX);
double gsl_blas_raw_dznrm2 (int N, const void * X, int incX);
double gsl_blas_raw_dzasum (int N, const void * X, int incX);

/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX gsl_blas_raw_isamax (int N, const float  X[], int incX);
CBLAS_INDEX gsl_blas_raw_idamax (int N, const double X[], int incX);
CBLAS_INDEX gsl_blas_raw_icamax (int N, const void * X, int incX);
CBLAS_INDEX gsl_blas_raw_izamax (int N, const void * X, int incX);


/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void gsl_blas_raw_sswap (int N,
                         float X[], int incX,
                         float Y[], int incY);

void gsl_blas_raw_scopy (int N,
                         const float X[], int incX,
                         float Y[], int incY);

void gsl_blas_raw_saxpy (int N, float alpha,
                         const float X[], int incX,
                         float Y[], int incY);

void gsl_blas_raw_dswap (int N,
                         double X[], int incX,
                         double Y[], int incY);

void gsl_blas_raw_dcopy (int N,
                         const double X[], int incX,
                         double Y[], int incY);

void gsl_blas_raw_daxpy (int N, double alpha,
                         const double X[], int incX, 
                         double Y[], int incY);

void gsl_blas_raw_cswap (int N,
                         void * X, int incX,
                         void * Y, int incY);

void gsl_blas_raw_ccopy (int N,
                         const void * X, int incX,
                         void * Y, int incY);

void gsl_blas_raw_caxpy (int N, const void * alpha,
                         const void * X, int incX,
                         void * Y, int incY);

void gsl_blas_raw_zswap (int N,
                         void * X, int incX,
                         void * Y, int incY);

void gsl_blas_raw_zcopy (int N,
                         const void * X, int incX,
                         void * Y, int incY);

void gsl_blas_raw_zaxpy (int N, const void * alpha,
                         const void * X, int incX,
                         void * Y, int incY);

/*
 * Routines with S and D prefix only
 */
void gsl_blas_raw_srotg (float a[], float b[], float c[], float s[]);

void gsl_blas_raw_srotmg (float d1[], float d2[], float b1[], float b2, float P[]);

void gsl_blas_raw_srot (int N,
                        float X[], int incX,
                        float Y[], int incY,
                        float c, float s);

void gsl_blas_raw_srotm (int N,
                         float X[], int incX,
                         float Y[], int incY,
                         const float P[]);

void gsl_blas_raw_drotg (double a[], double b[], double c[], double s[]);

void gsl_blas_raw_drotmg (double d1[], double d2[], double b1[],
                          double b2, double P[]);

void gsl_blas_raw_drot (int N,
                        double X[], int incX,
                        double Y[], int incY,
                        const double c, const double s);

void gsl_blas_raw_drotm (int N,
                         double X[], int incX,
                         double Y[], int incY,
			 const double P[]);

/*
 * Routines with S D C Z CS and ZD prefixes
 */
void gsl_blas_raw_sscal  (int N, float  alpha, float  X[], int incX);
void gsl_blas_raw_dscal  (int N, double alpha, double X[], int incX);
void gsl_blas_raw_cscal  (int N, const void * alpha, void * X, int incX);
void gsl_blas_raw_zscal  (int N, const void * alpha, void * X, int incX);
void gsl_blas_raw_csscal (int N, const float  alpha, void * X, int incX);
void gsl_blas_raw_zdscal (int N, const double alpha, void * X, int incX);


/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void gsl_blas_raw_sgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                         int M, int N,
                         const float alpha,
                         const float A[], int lda,
                         const float X[], int incX,
                         const float beta,
                         float Y[], int incY);

void gsl_blas_raw_sgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                         int M, int N,
                         int KL, int KU, const float alpha,
                         const float A[], int lda,
                         const float X[], int incX,
                         const float beta,
                         float Y[], int incY);

void gsl_blas_raw_strmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
                         const float A[], int lda,
                         float X[], int incX);

void gsl_blas_raw_stbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const float A[], int lda,
                         float X[], int incX);

void gsl_blas_raw_stpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const float Ap[],
			 float X[], int incX);

void gsl_blas_raw_strsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const float A[], int lda,
			 float X[], int incX);

void gsl_blas_raw_stbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const float A[], int lda,
                         float X[], int incX);

void gsl_blas_raw_stpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const float Ap[],
			 float X[], int incX);

void gsl_blas_raw_dgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
			 int M, int N,
                         const double alpha,
			 const double A[], int lda,
                         const double X[], int incX,
			 const double beta,
                         double Y[], int incY);

void gsl_blas_raw_dgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
			 int M, int N, int KL, int KU,
			 const double alpha,
                         const double A[], int lda,
			 const double X[], int incX,
			 const double beta,
			 double Y[], int incY);

void gsl_blas_raw_dtrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const double A[], int lda,
                         double X[], int incX);

void gsl_blas_raw_dtbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const double A[], int lda,
                         double X[], int incX);

void gsl_blas_raw_dtpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const double Ap[],
			 double X[], int incX);

void gsl_blas_raw_dtrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const double A[], int lda,
			 double X[], int incX);

void gsl_blas_raw_dtbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const double A[], int lda,
                         double X[], int incX);

void gsl_blas_raw_dtpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const double Ap[],
			 double X[], int incX);

void gsl_blas_raw_cgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         const void * X, int incX,
			 const void * beta,
                         void * Y, int incY);

void gsl_blas_raw_cgbmv (CBLAS_ORDER order,
                         CBLAS_TRANSPOSE TransA,
			 int M, int N, int KL, int KU,
			 const void * alpha,
                         const void * A, int lda,
			 const void * X, int incX,
			 const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_ctrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * A, int lda,
                         void * X, int incX);

void gsl_blas_raw_ctbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const void * A, int lda,
                         void * X, int incX);

void gsl_blas_raw_ctpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * Ap,
			 void * X, int incX);

void gsl_blas_raw_ctrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * A, int lda,
			 void * X, int incX);

void gsl_blas_raw_ctbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const void * A, int lda,
                         void * X, int incX);

void gsl_blas_raw_ctpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * Ap,
			 void * X, int incX);

void gsl_blas_raw_zgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         const void * X, int incX,
			 const void * beta,
                         void * Y, int incY);

void gsl_blas_raw_zgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
			 int M, int N, int KL, int KU,
			 const void * alpha,
                         const void * A, int lda,
			 const void * X, int incX,
			 const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_ztrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * A, int lda,
                         void * X, int incX);

void gsl_blas_raw_ztbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const void * A, int lda,
                         void * X, int incX);

void gsl_blas_raw_ztpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * Ap,
			 void * X, int incX);

void gsl_blas_raw_ztrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * A, int lda,
			 void * X, int incX);

void gsl_blas_raw_ztbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N, int K,
			 const void * A, int lda,
                         void * X, int incX);

void gsl_blas_raw_ztpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         int N,
			 const void * Ap,
			 void * X, int incX);

/*
 * Routines with S and D prefixes only
 */
void gsl_blas_raw_ssymv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const float alpha,
			 const float A[], int lda,
			 const float X[], int incX,
                         const float beta,
			 float Y[], int incY);

void gsl_blas_raw_ssbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N, int K,
			 const float alpha,
			 const float A[], int lda, 
			 const float X[], int incX,
                         const float beta,
			 float Y[], int incY);

void gsl_blas_raw_sspmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const float alpha,
			 const float Ap[],
                         const float X[], int incX,
                         const float beta,
			 float Y[], int incY);

void gsl_blas_raw_sger (CBLAS_ORDER order,
                        int M, int N,
                        const float alpha,
			const float X[], int incX,
                        const float Y[], int incY,
			float A[], int lda);

void gsl_blas_raw_ssyr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const float alpha,
			const float X[], int incX,
			float A[], int lda);

void gsl_blas_raw_sspr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const float alpha,
			const float X[], int incX,
			float Ap[]);

void gsl_blas_raw_ssyr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const float alpha,
			 const float X[], int incX, 
			 const float Y[], int incY,
			 float A[], int lda);

void gsl_blas_raw_sspr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const float alpha,
			 const float X[], int incX,
			 const float Y[], int incY,
			 float A[]);

void gsl_blas_raw_dsymv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const double alpha,
			 const double A[], int lda,
			 const double X[], int incX,
                         const double beta,
			 double Y[], int incY);

void gsl_blas_raw_dsbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N, int K,
			 const double alpha,
			 const double A[], int lda,
			 const double X[], int incX,
                         const double beta,
			 double Y[], int incY);

void gsl_blas_raw_dspmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const double alpha,
			 const double Ap[],
                         const double X[], int incX,
                         const double beta,
			 double Y[], int incY);

void gsl_blas_raw_dger (CBLAS_ORDER order,
                        int M, int N,
                        const double alpha,
			const double X[], int incX,
                        const double Y[], int incY,
			double A[], int lda);

void gsl_blas_raw_dsyr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const double alpha,
			const double X[], int incX,
			double A[], int lda);

void gsl_blas_raw_dspr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const double alpha,
			const double X[], int incX,
			double Ap[]);

void gsl_blas_raw_dsyr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const double alpha,
			 const double X[], int incX,
			 const double Y[], int incY,
			 double A[], int lda);

void gsl_blas_raw_dspr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const double alpha,
			 const double X[], int incX,
			 const double Y[], int incY,
			 double A[]);

/*
 * Routines with C and Z prefixes only
 */
void gsl_blas_raw_chemv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const void * alpha,
			 const void * A, int lda,
			 const void * X, int incX,
                         const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_chbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N, int K,
			 const void * alpha,
			 const void * A, int lda,
			 const void * X, int incX,
                         const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_chpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const void * alpha, const void * Ap,
                         const void * X, int incX,
                         const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_cgeru (CBLAS_ORDER order,
                         int M, int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * A, int lda);

void gsl_blas_raw_cgerc (CBLAS_ORDER order,
                         int M, int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * A, int lda);

void gsl_blas_raw_cher (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const float alpha,
			const void * X, int incX,
                        void * A, int lda);

void gsl_blas_raw_chpr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const float alpha[],
			const void * X, int incX,
			void * A);

void gsl_blas_raw_cher2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * A, int lda);

void gsl_blas_raw_chpr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * Ap);

void gsl_blas_raw_zhemv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const void * alpha,
			 const void * A, int lda,
			 const void * X, int incX,
                         const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_zhbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N, int K,
			 const void * alpha,
			 const void * A, int lda,
			 const void * X, int incX,
                         const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_zhpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
			 const void * alpha,
			 const void * Ap,
                         const void * X, int incX,
                         const void * beta,
			 void * Y, int incY);

void gsl_blas_raw_zgeru (CBLAS_ORDER order,
                         int M, int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * A, int lda);

void gsl_blas_raw_zgerc (CBLAS_ORDER order,
                         int M, int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * A, int lda);

void gsl_blas_raw_zher (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const double alpha,
			const void * X, int incX,
                        void * A, int lda);

void gsl_blas_raw_zhpr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                        int N,
			const double alpha[],
			const void * X, int incX,
			void * A);

void gsl_blas_raw_zher2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * A, int lda);

void gsl_blas_raw_zhpr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                         int N,
                         const void * alpha,
			 const void * X, int incX,
                         const void * Y, int incY,
			 void * Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void gsl_blas_raw_sgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 int M, int N, int K,
			 const float alpha,
			 const float A[], int lda,
			 const float B[], int ldb,
                         const float beta,
			 float C[], int ldc);

void gsl_blas_raw_ssymm (CBLAS_ORDER Order, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
			 int M, int N,
                         const float alpha,
			 const float A[], int lda,
                         const float B[], int ldb,
			 const float beta,
                         float C[], int ldc);

void gsl_blas_raw_ssyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
			 int N, int K,
                         const float alpha,
			 const float A[], int lda,
                         const float beta,
			 float C[], int ldc);

void gsl_blas_raw_ssyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
			  int N, int K,
                          const float alpha,
			  const float A[], int lda,
                          const float B[], int ldb,
			  const float beta,
                          float C[], int ldc);

void gsl_blas_raw_strmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const float alpha,
			 const float A[], int lda,
                         float B[], int ldb);

void gsl_blas_raw_strsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const float alpha,
			 const float A[], int lda,
                         float B[], int ldb);

void gsl_blas_raw_dgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 int M, int N, int K,
			 const double alpha,
			 const double A[], int lda,
			 const double B[], int ldb,
                         const double beta,
			 double C[], int ldc);

void gsl_blas_raw_dsymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 int M, int N,
                         const double alpha,
			 const double A[], int lda,
                         const double B[], int ldb,
			 const double beta,
                         double C[], int ldc);

void gsl_blas_raw_dsyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 int N, int K,
                         const double alpha,
			 const double A[], int lda,
                         const double beta,
			 double C[], int ldc);

void gsl_blas_raw_dsyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  int N, int K,
                          const double alpha,
			  const double A[], int lda,
                          const double B[], int ldb,
			  const double beta,
                          double C[], int ldc);

void gsl_blas_raw_dtrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const double alpha,
			 const double A[], int lda,
                         double B[], int ldb);

void gsl_blas_raw_dtrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const double alpha,
			 const double A[], int lda,
                         double B[], int ldb);

void gsl_blas_raw_cgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 int M, int N, int K,
			 const void * alpha,
			 const void * A, int lda, 
			 const void * B, int ldb,
                         const void * beta,
			 void * C, int ldc);

void gsl_blas_raw_csymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         const void * B, int ldb,
			 const void * beta,
                         void * C, int ldc);

void gsl_blas_raw_csyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 int N, int K,
                         const void * alpha,
			 const void * A, int lda,
                         const void * beta,
			 void * C, int ldc);

void gsl_blas_raw_csyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  int N, int K,
                          const void * alpha,
			  const void * A, int lda,
                          const void * B, int ldb,
			  const void * beta,
                          void * C, int ldc);

void gsl_blas_raw_ctrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         void * B, int ldb);

void gsl_blas_raw_ctrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         void * B, int ldb);

void gsl_blas_raw_zgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 int M, int N, int K,
			 const void * alpha,
			 const void * A, int lda,
			 const void * B, int ldb,
                         const void * beta,
			 void * C, int ldc);

void gsl_blas_raw_zsymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         const void * B, int ldb,
			 const void * beta,
                         void * C, int ldc);

void gsl_blas_raw_zsyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 int N, int K,
                         const void * alpha,
			 const void * A, int lda,
                         const void * beta,
			 void * C, int ldc);

void gsl_blas_raw_zsyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  int N, int K,
                          const void * alpha,
			  const void * A, int lda,
                          const void * B, int ldb,
			  const void * beta,
                          void * C, int ldc);

void gsl_blas_raw_ztrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         void * B, int ldb);

void gsl_blas_raw_ztrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         void * B, int ldb);

/*
 * Routines with prefixes C and Z only
 */
void gsl_blas_raw_chemm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         const void * B, int ldb,
			 const void * beta,
                         void * C, int ldc);

void gsl_blas_raw_cherk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 int N, int K,
                         const float alpha,
			 const void * A, int lda,
                         const float beta,
			 void * C, int ldc);

void gsl_blas_raw_cher2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  int N, int K,
                          const void * alpha,
			  const void * A, int lda,
                          const void * B, int ldb,
			  const float beta,
                          void * C, int ldc);

void gsl_blas_raw_zhemm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 int M, int N,
                         const void * alpha,
			 const void * A, int lda,
                         const void * B, int ldb,
			 const void * beta,
                         void * C, int ldc);

void gsl_blas_raw_zherk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 int N, int K,
                         const double alpha,
			 const void * A, int lda,
                         const double beta,
			 void * C, int ldc);

void gsl_blas_raw_zher2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  int N, int K,
                          const void * alpha,
			  const void * A, int lda,
                          const void * B, int ldb,
			  const double beta,
                          void * C, int ldc);


#if defined(HAVE_INLINE) && defined(HAVE_CBLAS)

/* insert inline cblas implementation of above here */

#endif /* defined(HAVE_INLINE) && defined(HAVE_CBLAS) */


#endif /* !GSL_BLAS_RAW_H_ */
