/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Prototypes for level 3 BLAS functions.
 * Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef GSL_BLAS_RAW_L3_H
#define GSL_BLAS_RAW_L3_H

#include <gsl_complex.h>
#include <gsl_blas_types.h>


/* GEMM */

void gsl_blas_raw_sgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 size_t M, size_t N, size_t K,
			 float alpha,
			 const float A[], int lda,
			 const float B[], int ldb,
                         float beta,
			 float C[], int ldc);

void gsl_blas_raw_dgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 size_t M, size_t N, size_t K,
			 double alpha,
			 const double A[], int lda,
			 const double B[], int ldb,
                         double beta,
			 double C[], int ldc);

void gsl_blas_raw_cgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 size_t M, size_t N, size_t K,
			 const gsl_complex_packed_float alpha,
			 const gsl_complex_packed_array_float A, int lda, 
			 const gsl_complex_packed_array_float B, int ldb,
                         const gsl_complex_packed_float beta,
			 gsl_complex_packed_array_float C, int ldc);

void gsl_blas_raw_zgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
			 size_t M, size_t N, size_t K,
			 const gsl_complex_packed alpha,
			 const gsl_complex_packed_array A, int lda,
			 const gsl_complex_packed_array B, int ldb,
                         const gsl_complex_packed beta,
			 gsl_complex_packed_array C, int ldc);


/* SYMM */

void gsl_blas_raw_ssymm (CBLAS_SIDE Side, CBLAS_UPLO Uplo,
			 size_t M, size_t N,
                         float alpha,
			 const float A[], int lda,
                         const float B[], int ldb,
			 float beta,
                         float C[], int ldc);

void gsl_blas_raw_dsymm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 size_t M, size_t N,
                         double alpha,
			 const double A[], int lda,
                         const double B[], int ldb,
			 double beta,
                         double C[], int ldc);

void gsl_blas_raw_csymm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
			 const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_array_float B, int ldb,
			 const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc);

void gsl_blas_raw_zsymm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 size_t M, size_t N,
                         const gsl_complex_packed alpha,
			 const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array B, int ldb,
			 const gsl_complex_packed beta,
                         gsl_complex_packed_array C, int ldc);


/* SYRK */

void gsl_blas_raw_ssyrk (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
			 size_t N, size_t K,
                         float alpha,
			 const float A[], int lda,
                         float beta,
			 float C[], int ldc);

void gsl_blas_raw_dsyrk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 size_t N, size_t K,
                         double alpha,
			 const double A[], int lda,
                         double beta,
			 double C[], int ldc);

void gsl_blas_raw_csyrk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 size_t N, size_t K,
                         const gsl_complex_packed_float alpha,
			 const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_float beta,
			 gsl_complex_packed_array_float C, int ldc);

void gsl_blas_raw_zsyrk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 size_t N, size_t K,
                         const gsl_complex_packed alpha,
			 const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed beta,
			 gsl_complex_packed_array C, int ldc);


/* SYR2K */

void gsl_blas_raw_ssyr2k (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
			  size_t N, size_t K,
                          float alpha,
			  const float A[], int lda,
                          const float B[], int ldb,
			  float beta,
                          float C[], int ldc);

void gsl_blas_raw_dsyr2k (CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  size_t N, size_t K,
                          double alpha,
			  const double A[], int lda,
                          const double B[], int ldb,
			  double beta,
                          double C[], int ldc);

void gsl_blas_raw_csyr2k (CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  size_t N, size_t K,
                          const gsl_complex_packed_float alpha,
			  const gsl_complex_packed_array_float A, int lda,
                          const gsl_complex_packed_array_float B, int ldb,
			  const gsl_complex_packed_float beta,
                          gsl_complex_packed_array_float C, int ldc);

void gsl_blas_raw_zsyr2k (CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  size_t N, size_t K,
                          const gsl_complex_packed alpha,
			  const gsl_complex_packed_array A, int lda,
                          const gsl_complex_packed_array B, int ldb,
			  const gsl_complex_packed beta,
                          gsl_complex_packed_array C, int ldc);


/* TRMM */

void gsl_blas_raw_strmm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         float alpha,
			 const float A[], int lda,
                         float B[], int ldb);

void gsl_blas_raw_dtrmm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         double alpha,
			 const double A[], int lda,
                         double B[], int ldb);

void gsl_blas_raw_ctrmm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
			 const gsl_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float B, int ldb);

void gsl_blas_raw_ztrmm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         const gsl_complex_packed alpha,
			 const gsl_complex_packed_array A, int lda,
                         gsl_complex_packed_array B, int ldb);


/* TRSM */

void gsl_blas_raw_strsm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         float alpha,
			 const float A[], int lda,
                         float B[], int ldb);

void gsl_blas_raw_dtrsm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         double alpha,
			 const double A[], int lda,
                         double B[], int ldb);

void gsl_blas_raw_ctrsm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
			 const gsl_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float B, int ldb);

void gsl_blas_raw_ztrsm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
			 size_t M, size_t N,
                         const gsl_complex_packed alpha,
			 const gsl_complex_packed_array A, int lda,
                         gsl_complex_packed_array B, int ldb);


/* HEMM */

void gsl_blas_raw_chemm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
			 const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_array_float B, int ldb,
			 const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc);

void gsl_blas_raw_zhemm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
			 size_t M, size_t N,
                         const gsl_complex_packed alpha,
			 const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array B, int ldb,
			 const gsl_complex_packed beta,
                         gsl_complex_packed_array C, int ldc);


/* HERK */

void gsl_blas_raw_cherk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 size_t N, size_t K,
                         float alpha,
			 const gsl_complex_packed_array_float A, int lda,
                         float beta,
			 gsl_complex_packed_array_float C, int ldc);

void gsl_blas_raw_zherk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
			 size_t N, size_t K,
                         double alpha,
			 const gsl_complex_packed_array A, int lda,
                         double beta,
			 gsl_complex_packed_array C, int ldc);


/* HER2K */

void gsl_blas_raw_cher2k (CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  size_t N, size_t K,
                          const gsl_complex_packed_float alpha,
			  const gsl_complex_packed_array_float A, int lda,
                          const gsl_complex_packed_array_float B, int ldb,
			  float beta,
                          gsl_complex_packed_array_float C, int ldc);


void gsl_blas_raw_zher2k (CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
			  size_t N, size_t K,
                          const gsl_complex_packed alpha,
			  const gsl_complex_packed_array A, int lda,
                          const gsl_complex_packed_array B, int ldb,
			  double beta,
                          gsl_complex_packed_array C, int ldc);


#endif /* !GSL_BLAS_RAW_L3_H */
