/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Prototypes for level 2 BLAS functions.
 * Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef GSL_BLAS_RAW_L2_H
#define GSL_BLAS_RAW_L2_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas_types.h>


/* GEMV */

void gsl_blas_raw_sgemv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N,
                         float alpha,
                         const float A[], int lda,
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY);

void gsl_blas_raw_dgemv (CBLAS_TRANSPOSE_t TransA,
			 size_t M, size_t N,
                         double alpha,
			 const double A[], int lda,
                         const double X[], size_t incX,
			 double beta,
                         double Y[], size_t incY);

void gsl_blas_raw_cgemv (CBLAS_TRANSPOSE_t TransA,
			 size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
			 const gsl_const_complex_packed_array_float A, int lda,
                         const gsl_const_complex_packed_array_float X, size_t incX,
			 const gsl_const_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY);

void gsl_blas_raw_zgemv (CBLAS_TRANSPOSE_t TransA,
			 size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array A, int lda,
                         const gsl_const_complex_packed_array X, size_t incX,
			 const gsl_const_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY);

/* GBMV */

void gsl_blas_raw_sgbmv (CBLAS_TRANSPOSE_t TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
			 float alpha,
                         const float A[], int lda,
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY);

void gsl_blas_raw_dgbmv (CBLAS_TRANSPOSE_t TransA,
			 size_t M, size_t N, size_t KL, size_t KU,
			 double alpha,
                         const double A[], int lda,
			 const double X[], size_t incX,
			 double beta,
			 double Y[], size_t incY);

void gsl_blas_raw_cgbmv (CBLAS_TRANSPOSE_t TransA,
			 size_t M, size_t N, size_t KL, size_t KU,
			 const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float A, int lda,
			 const gsl_const_complex_packed_array_float X, size_t incX,
			 const gsl_const_complex_packed_float beta,
			 gsl_complex_packed_array_float Y, size_t incY);

void gsl_blas_raw_zgbmv (CBLAS_TRANSPOSE_t TransA,
			 size_t M, size_t N, size_t KL, size_t KU,
			 const gsl_const_complex_packed alpha,
                         const gsl_const_complex_packed_array A, int lda,
			 const gsl_const_complex_packed_array X, size_t incX,
			 const gsl_const_complex_packed beta,
			 gsl_complex_packed_array Y, size_t incY);


/* TRMV */

void gsl_blas_raw_strmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
                         const float A[], int lda,
                         float X[], size_t incX);

void gsl_blas_raw_dtrmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const double A[], int lda,
                         double X[], size_t incX);

void gsl_blas_raw_ctrmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX);

void gsl_blas_raw_ztrmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX);


/* TBMV */

void gsl_blas_raw_stbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const float A[], int lda,
                         float X[], size_t incX);

void gsl_blas_raw_dtbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const double A[], int lda,
                         double X[], size_t incX);

void gsl_blas_raw_ctbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX);

void gsl_blas_raw_ztbmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX);


/* TPMV */

void gsl_blas_raw_stpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const float Ap[],
			 float X[], size_t incX);

void gsl_blas_raw_dtpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const double Ap[],
			 double X[], size_t incX);

void gsl_blas_raw_ctpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array_float Ap,
			 gsl_complex_packed_array_float X, size_t incX);

void gsl_blas_raw_ztpmv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array Ap,
			 gsl_complex_packed_array X, size_t incX);

/* TRSV */

void gsl_blas_raw_strsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const float A[], int lda,
			 float X[], size_t incX);

void gsl_blas_raw_dtrsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const double A[], int lda,
			 double X[], size_t incX);

void gsl_blas_raw_ctrsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array_float A, int lda,
			 gsl_complex_packed_array_float X, size_t incX);

void gsl_blas_raw_ztrsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array A, int lda,
			 gsl_complex_packed_array X, size_t incX);


/* TBSV */

void gsl_blas_raw_stbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const float A[], int lda,
                         float X[], size_t incX);

void gsl_blas_raw_dtbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const double A[], int lda,
                         double X[], size_t incX);

void gsl_blas_raw_ctbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const gsl_const_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX);

void gsl_blas_raw_ztbsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N, size_t K,
			 const gsl_const_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX);


/* TPSV */

void gsl_blas_raw_stpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const float Ap[],
			 float X[], size_t incX);

void gsl_blas_raw_dtpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const double Ap[],
			 double X[], size_t incX);

void gsl_blas_raw_ctpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array_float Ap,
			 gsl_complex_packed_array_float X, size_t incX);

void gsl_blas_raw_ztpsv (CBLAS_UPLO_t Uplo,
                         CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         size_t N,
			 const gsl_const_complex_packed_array Ap,
			 gsl_complex_packed_array X, size_t incX);


/* SYMV */

void gsl_blas_raw_ssymv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 float alpha,
			 const float A[], int lda,
			 const float X[], size_t incX,
                         float beta,
			 float Y[], size_t incY);

void gsl_blas_raw_dsymv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 double alpha,
			 const double A[], int lda,
			 const double X[], size_t incX,
                         double beta,
			 double Y[], size_t incY);


/* SBMV */

void gsl_blas_raw_ssbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
			 float alpha,
			 const float A[], int lda, 
			 const float X[], size_t incX,
                         float beta,
			 float Y[], size_t incY);

void gsl_blas_raw_dsbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
			 double alpha,
			 const double A[], int lda,
			 const double X[], size_t incX,
                         double beta,
			 double Y[], size_t incY);

/* SPMV */

void gsl_blas_raw_sspmv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 float alpha,
			 const float Ap[],
                         const float X[], size_t incX,
                         float beta,
			 float Y[], size_t incY);

void gsl_blas_raw_dspmv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 double alpha,
			 const double Ap[],
                         const double X[], size_t incX,
                         double beta,
			 double Y[], size_t incY);

/* GER */

void gsl_blas_raw_sger (size_t M, size_t N,
                        float alpha,
			const float X[], size_t incX,
                        const float Y[], size_t incY,
			float A[], int lda);

void gsl_blas_raw_dger (size_t M, size_t N,
                        double alpha,
			const double X[], size_t incX,
                        const double Y[], size_t incY,
			double A[], int lda);


/* SYR */

void gsl_blas_raw_ssyr (CBLAS_UPLO_t Uplo,
                        size_t N,
			float alpha,
			const float X[], size_t incX,
			float A[], int lda);

void gsl_blas_raw_dsyr (CBLAS_UPLO_t Uplo,
                        size_t N,
			double alpha,
			const double X[], size_t incX,
			double A[], int lda);


/* SPR */

void gsl_blas_raw_sspr (CBLAS_UPLO_t Uplo,
                        size_t N,
			float alpha,
			const float X[], size_t incX,
			float Ap[]);

void gsl_blas_raw_dspr (CBLAS_UPLO_t Uplo,
                        size_t N,
			double alpha,
			const double X[], size_t incX,
			double Ap[]);


/* SYR2 */

void gsl_blas_raw_ssyr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
			 float alpha,
			 const float X[], size_t incX, 
			 const float Y[], size_t incY,
			 float A[], int lda);

void gsl_blas_raw_dsyr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
			 double alpha,
			 const double X[], size_t incX,
			 const double Y[], size_t incY,
			 double A[], int lda);


/* SPR2 */

void gsl_blas_raw_sspr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
			 float alpha,
			 const float X[], size_t incX,
			 const float Y[], size_t incY,
			 float A[]);

void gsl_blas_raw_dspr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
			 double alpha,
			 const double X[], size_t incX,
			 const double Y[], size_t incY,
			 double A[]);


/* HEMV */

void gsl_blas_raw_chemv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 const gsl_const_complex_packed_float alpha,
			 const gsl_const_complex_packed_array_float A, int lda,
			 const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
			 gsl_complex_packed_array_float Y, size_t incY);

void gsl_blas_raw_zhemv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array A, int lda,
			 const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
			 gsl_complex_packed_array Y, size_t incY);


/* HBMV */

void gsl_blas_raw_chbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
			 const gsl_const_complex_packed_float alpha,
			 const gsl_const_complex_packed_array_float A, int lda,
			 const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
			 gsl_complex_packed_array_float Y, size_t incY);

void gsl_blas_raw_zhbmv (CBLAS_UPLO_t Uplo,
                         size_t N, size_t K,
			 const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array A, int lda,
			 const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
			 gsl_complex_packed_array Y, size_t incY);


/* HPMV */

void gsl_blas_raw_chpmv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 const gsl_const_complex_packed_float alpha,
                         const gsl_const_complex_packed_array_float Ap,
                         const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_float beta,
			 gsl_complex_packed_array_float Y, size_t incY);

void gsl_blas_raw_zhpmv (CBLAS_UPLO_t Uplo,
                         size_t N,
			 const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array Ap,
                         const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed beta,
			 gsl_complex_packed_array Y, size_t incY);


/* GERU */

void gsl_blas_raw_cgeru (size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
			 const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
			 gsl_complex_packed_array_float A, int lda);

void gsl_blas_raw_zgeru (size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
			 gsl_complex_packed_array A, int lda);


/* GERC */

void gsl_blas_raw_cgerc (size_t M, size_t N,
                         const gsl_const_complex_packed_float alpha,
			 const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
			 gsl_complex_packed_array_float A, int lda);

void gsl_blas_raw_zgerc (size_t M, size_t N,
                         const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
			 gsl_complex_packed_array A, int lda);

/* HER */

void gsl_blas_raw_cher (CBLAS_UPLO_t Uplo,
                        size_t N,
			float alpha,
			const gsl_const_complex_packed_array_float X, size_t incX,
                        gsl_complex_packed_array_float A, int lda);

void gsl_blas_raw_zher (CBLAS_UPLO_t Uplo,
                        size_t N,
			double alpha,
			const gsl_const_complex_packed_array X, size_t incX,
                        gsl_complex_packed_array A, int lda);


/* HPR */

void gsl_blas_raw_chpr (CBLAS_UPLO_t Uplo,
                        size_t N,
			float alpha,
			const gsl_const_complex_packed_array_float X, size_t incX,
			gsl_complex_packed_array_float A);

void gsl_blas_raw_zhpr (CBLAS_UPLO_t Uplo,
                        size_t N,
			double alpha,
			const gsl_const_complex_packed_array X, size_t incX,
			gsl_complex_packed_array A);


/* HER2 */

void gsl_blas_raw_cher2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed_float alpha,
			 const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
			 gsl_complex_packed_array_float A, int lda);

void gsl_blas_raw_zher2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
			 gsl_complex_packed_array A, int lda);


/* HPR2 */

void gsl_blas_raw_chpr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed_float alpha,
			 const gsl_const_complex_packed_array_float X, size_t incX,
                         const gsl_const_complex_packed_array_float Y, size_t incY,
			 gsl_complex_packed_array_float Ap);

void gsl_blas_raw_zhpr2 (CBLAS_UPLO_t Uplo,
                         size_t N,
                         const gsl_const_complex_packed alpha,
			 const gsl_const_complex_packed_array X, size_t incX,
                         const gsl_const_complex_packed_array Y, size_t incY,
			 gsl_complex_packed_array Ap);


#endif /* !GSL_BLAS_RAW_L2_H */
