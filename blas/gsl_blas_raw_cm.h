/* blas/gsl_blas_raw_cm.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* Raw BLAS interface for column-major matrices.
 * Based on draft BLAST C interface specification  [Jul 7 1998]
 */
#ifndef __GSL_BLAS_RAW_CM_H__
#define __GSL_BLAS_RAW_CM_H__

#include <gsl/gsl_blas_raw_L1.h>
#include <gsl/gsl_blas_types.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */


/* GEMV */
void gsl_blas_raw_sgemv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_raw_dgemv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

void gsl_blas_raw_cgemv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_raw_zgemv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

/* GBMV */

void gsl_blas_raw_sgbmv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_raw_dgbmv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

void gsl_blas_raw_cgbmv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_raw_zgbmv_cm (CBLAS_TRANSPOSE_t TransA,
                            size_t M, size_t N, size_t KL, size_t KU,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* TRMV */

void gsl_blas_raw_strmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_raw_dtrmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_raw_ctrmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_raw_ztrmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);


/* TBMV */

void gsl_blas_raw_stbmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_raw_dtbmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_raw_ctbmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_raw_ztbmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);


/* TPMV */

void gsl_blas_raw_stpmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const float Ap[],
                            float X[], int incX);

void gsl_blas_raw_dtpmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const double Ap[],
                            double X[], int incX);

void gsl_blas_raw_ctpmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);

void gsl_blas_raw_ztpmv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);

/* TRSV */

void gsl_blas_raw_strsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_raw_dtrsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_raw_ctrsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_raw_ztrsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * A, int lda,
                            void * X, int incX);


/* TBSV */

void gsl_blas_raw_stbsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const float A[], int lda,
                            float X[], int incX);

void gsl_blas_raw_dtbsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const double A[], int lda,
                            double X[], int incX);

void gsl_blas_raw_ctbsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);

void gsl_blas_raw_ztbsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N, size_t K,
                            const void * A, int lda,
                            void * X, int incX);


/* TPSV */

void gsl_blas_raw_stpsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const float Ap[],
                            float X[], int incX);

void gsl_blas_raw_dtpsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const double Ap[],
                            double X[], int incX);

void gsl_blas_raw_ctpsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);

void gsl_blas_raw_ztpsv_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                            size_t N,
                            const void * Ap,
                            void * X, int incX);


/* SYMV */

void gsl_blas_raw_ssymv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_raw_dsymv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);


/* SBMV */

void gsl_blas_raw_ssbmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N, size_t K,
                            float alpha,
                            const float A[], int lda, 
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_raw_dsbmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

/* SPMV */

void gsl_blas_raw_sspmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            float alpha,
                            const float Ap[],
                            const float X[], int incX,
                            float beta,
                            float Y[], int incY);

void gsl_blas_raw_dspmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            double alpha,
                            const double Ap[],
                            const double X[], int incX,
                            double beta,
                            double Y[], int incY);

/* GER */

void gsl_blas_raw_sger_cm (size_t M, size_t N,
                           float alpha,
                           const float X[], int incX,
                           const float Y[], int incY,
                           float A[], int lda);

void gsl_blas_raw_dger_cm (size_t M, size_t N,
                           double alpha,
                           const double X[], int incX,
                           const double Y[], int incY,
                           double A[], int lda);


/* SYR */

void gsl_blas_raw_ssyr_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           float alpha,
                           const float X[], int incX,
                           float A[], int lda);

void gsl_blas_raw_dsyr_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           double alpha,
                           const double X[], int incX,
                           double A[], int lda);


/* SPR */

void gsl_blas_raw_sspr_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           float alpha,
                           const float X[], int incX,
                           float Ap[]);

void gsl_blas_raw_dspr_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           double alpha,
                           const double X[], int incX,
                           double Ap[]);


/* SYR2 */

void gsl_blas_raw_ssyr2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            float alpha,
                            const float X[], int incX, 
                            const float Y[], int incY,
                            float A[], int lda);

void gsl_blas_raw_dsyr2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            double alpha,
                            const double X[], int incX,
                            const double Y[], int incY,
                            double A[], int lda);


/* SPR2 */

void gsl_blas_raw_sspr2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            float alpha,
                            const float X[], int incX,
                            const float Y[], int incY,
                            float A[]);

void gsl_blas_raw_dspr2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            double alpha,
                            const double X[], int incX,
                            const double Y[], int incY,
                            double A[]);


/* HEMV */

void gsl_blas_raw_chemv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_raw_zhemv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* HBMV */

void gsl_blas_raw_chbmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_raw_zhbmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* HPMV */

void gsl_blas_raw_chpmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha, const void * Ap,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);

void gsl_blas_raw_zhpmv_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha,
                            const void * Ap,
                            const void * X, int incX,
                            const void * beta,
                            void * Y, int incY);


/* GERU */

void gsl_blas_raw_cgeru_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

void gsl_blas_raw_zgeru_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);


/* GERC */

void gsl_blas_raw_cgerc_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

void gsl_blas_raw_zgerc_cm (size_t M, size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

/* HER */

void gsl_blas_raw_cher_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           float alpha,
                           const void * X, int incX,
                           void * A, int lda);

void gsl_blas_raw_zher_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           double alpha,
                           const void * X, int incX,
                           void * A, int lda);


/* HPR */

void gsl_blas_raw_chpr_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           float alpha,
                           const void * X, int incX,
                           void * A);

void gsl_blas_raw_zhpr_cm (CBLAS_UPLO_t Uplo,
                           size_t N,
                           double alpha,
                           const void * X, int incX,
                           void * A);


/* HER2 */

void gsl_blas_raw_cher2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);

void gsl_blas_raw_zher2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * A, int lda);


/* HPR2 */

void gsl_blas_raw_chpr2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * Ap);

void gsl_blas_raw_zhpr2_cm (CBLAS_UPLO_t Uplo,
                            size_t N,
                            const void * alpha,
                            const void * X, int incX,
                            const void * Y, int incY,
                            void * Ap);


/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/* GEMM */

void gsl_blas_raw_sgemm_cm (CBLAS_TRANSPOSE_t TransA,
                            CBLAS_TRANSPOSE_t TransB,
                            size_t M, size_t N, size_t K,
                            float alpha,
                            const float A[], int lda,
                            const float B[], int ldb,
                            float beta,
                            float C[], int ldc);

void gsl_blas_raw_dgemm_cm (CBLAS_TRANSPOSE_t TransA,
                            CBLAS_TRANSPOSE_t TransB,
                            size_t M, size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            const double B[], int ldb,
                            double beta,
                            double C[], int ldc);

void gsl_blas_raw_cgemm_cm (CBLAS_TRANSPOSE_t TransA,
                            CBLAS_TRANSPOSE_t TransB,
                            size_t M, size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda, 
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_raw_zgemm_cm (CBLAS_TRANSPOSE_t TransA,
                            CBLAS_TRANSPOSE_t TransB,
                            size_t M, size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);


/* SYMM */

void gsl_blas_raw_ssymm_cm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            const float B[], int ldb,
                            float beta,
                            float C[], int ldc);

void gsl_blas_raw_dsymm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            const double B[], int ldb,
                            double beta,
                            double C[], int ldc);

void gsl_blas_raw_csymm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_raw_zsymm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);


/* SYRK */

void gsl_blas_raw_ssyrk_cm (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                            size_t N, size_t K,
                            float alpha,
                            const float A[], int lda,
                            float beta,
                            float C[], int ldc);

void gsl_blas_raw_dsyrk_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t Trans,
                            size_t N, size_t K,
                            double alpha,
                            const double A[], int lda,
                            double beta,
                            double C[], int ldc);

void gsl_blas_raw_csyrk_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t Trans,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_raw_zsyrk_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t Trans,
                            size_t N, size_t K,
                            const void * alpha,
                            const void * A, int lda,
                            const void * beta,
                            void * C, int ldc);


/* SYR2K */

void gsl_blas_raw_ssyr2k_cm (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                             size_t N, size_t K,
                             float alpha,
                             const float A[], int lda,
                             const float B[], int ldb,
                             float beta,
                             float C[], int ldc);

void gsl_blas_raw_dsyr2k_cm (CBLAS_UPLO_t Uplo,
                             CBLAS_TRANSPOSE_t Trans,
                             size_t N, size_t K,
                             double alpha,
                             const double A[], int lda,
                             const double B[], int ldb,
                             double beta,
                             double C[], int ldc);

void gsl_blas_raw_csyr2k_cm (CBLAS_UPLO_t Uplo,
                             CBLAS_TRANSPOSE_t Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             const void * beta,
                             void * C, int ldc);

void gsl_blas_raw_zsyr2k_cm (CBLAS_UPLO_t Uplo,
                             CBLAS_TRANSPOSE_t Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             const void * beta,
                             void * C, int ldc);


/* TRMM */

void gsl_blas_raw_strmm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            float B[], int ldb);

void gsl_blas_raw_dtrmm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            double B[], int ldb);

void gsl_blas_raw_ctrmm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);

void gsl_blas_raw_ztrmm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);


/* TRSM */

void gsl_blas_raw_strsm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            float alpha,
                            const float A[], int lda,
                            float B[], int ldb);

void gsl_blas_raw_dtrsm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            double alpha,
                            const double A[], int lda,
                            double B[], int ldb);

void gsl_blas_raw_ctrsm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);

void gsl_blas_raw_ztrsm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                            CBLAS_DIAG_t Diag,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            void * B, int ldb);


/* HEMM */

void gsl_blas_raw_chemm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);

void gsl_blas_raw_zhemm_cm (CBLAS_SIDE_t Side,
                            CBLAS_UPLO_t Uplo,
                            size_t M, size_t N,
                            const void * alpha,
                            const void * A, int lda,
                            const void * B, int ldb,
                            const void * beta,
                            void * C, int ldc);


/* HERK */

void gsl_blas_raw_cherk_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t Trans,
                            size_t N, size_t K,
                            float alpha,
                            const void * A, int lda,
                            float beta,
                            void * C, int ldc);

void gsl_blas_raw_zherk_cm (CBLAS_UPLO_t Uplo,
                            CBLAS_TRANSPOSE_t Trans,
                            size_t N, size_t K,
                            double alpha,
                            const void * A, int lda,
                            double beta,
                            void * C, int ldc);


/* HER2K */

void gsl_blas_raw_cher2k_cm (CBLAS_UPLO_t Uplo,
                             CBLAS_TRANSPOSE_t Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             float beta,
                             void * C, int ldc);


void gsl_blas_raw_zher2k_cm (CBLAS_UPLO_t Uplo,
                             CBLAS_TRANSPOSE_t Trans,
                             size_t N, size_t K,
                             const void * alpha,
                             const void * A, int lda,
                             const void * B, int ldb,
                             double beta,
                             void * C, int ldc);


#if defined(HAVE_INLINE) && defined(HAVE_CBLAS)
#include <cblas.h>

/* insert inline cblas implementation of above here */

#endif /* defined(HAVE_INLINE) && defined(HAVE_CBLAS) */


__END_DECLS

#endif /* __GSL_BLAS_RAW_CM_H__ */
