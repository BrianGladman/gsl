/*
 * Author:  G. Jungman  and  B. Gough
 * RCS:     $Id$
 */
/* Native implementation of row-major operations.
 * Conforms to gsl_blas_raw interface.
 *
 */
#include <gsl_math.h>
#include "complex_internal.h"
#include "gsl_blas_raw.h"



/* ===========================================================================
 * level 1 BLAS functions
 * ===========================================================================
 */

float  gsl_blas_raw_sdsdot (size_t N,
                            float alpha,
                            const float X[], size_t incX,
                            const float Y[], size_t incY)
{
#define INIT_VAL  alpha
#define ACC_TYPE  float
#define BASE_TYPE float
#include "source_dot_r.h"
#undef ACC_TYPE 
#undef BASE_TYPE
#undef INIT_VAL
}

double gsl_blas_raw_dsdot (size_t N,
                           const float X[], size_t incX,
                           const float Y[], size_t incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  double
#define BASE_TYPE float
#include "source_dot_r.h"
#undef ACC_TYPE 
#undef BASE_TYPE
#undef INIT_VAL
}

float  gsl_blas_raw_sdot (size_t N,
                          const float X[], size_t incX,
                          const float Y[], size_t incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  float
#define BASE_TYPE float
#include "source_dot_r.h"
#undef ACC_TYPE 
#undef BASE_TYPE
#undef INIT_VAL
}
                          
double gsl_blas_raw_ddot (size_t N,
                          const double X[], size_t incX,
                          const double Y[], size_t incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  double
#define BASE_TYPE double
#include "source_dot_r.h"
#undef ACC_TYPE 
#undef BASE_TYPE
#undef INIT_VAL
}


void gsl_blas_raw_cdotu (size_t N,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_float result)
{
#define BASE_TYPE float
#define CONJ_SIGN 1.0
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE_TYPE
}

void gsl_blas_raw_cdotc (size_t N,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_float result)
{
#define BASE_TYPE float
#define CONJ_SIGN (-1.0)
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE_TYPE
}

void gsl_blas_raw_zdotu (size_t N,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed_array Y, size_t incY,
                         gsl_complex_packed result)
{
#define BASE_TYPE double
#define CONJ_SIGN 1.0
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE_TYPE
}

void gsl_blas_raw_zdotc (size_t N,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed_array Y, size_t incY,
                         gsl_complex_packed result)
{
#define BASE_TYPE double
#define CONJ_SIGN (-1.0)
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE_TYPE
}


float  gsl_blas_raw_snrm2  (size_t N, const float  X[], size_t incX)
{
#define BASE_TYPE float
#include "source_nrm2_r.h"
#undef BASE_TYPE
}

double gsl_blas_raw_dnrm2  (size_t N, const double X[], size_t incX)
{
#define BASE_TYPE double
#include "source_nrm2_r.h"
#undef BASE_TYPE
}

float  gsl_blas_raw_scnrm2 (size_t N, const gsl_complex_packed_array_float X, size_t incX)
{
#define BASE_TYPE float
#include "source_nrm2_c.h"
#undef BASE_TYPE
}

double gsl_blas_raw_dznrm2 (size_t N, const gsl_complex_packed_array X, size_t incX)
{
#define BASE_TYPE double
#include "source_nrm2_c.h"
#undef BASE_TYPE
}

float  gsl_blas_raw_sasum (size_t N, const float X[], size_t incX)
{
#define BASE_TYPE float
#include "source_asum_r.h"
#undef BASE_TYPE
}

double gsl_blas_raw_dasum (size_t N, const double X[], size_t incX)
{
#define BASE_TYPE double
#include "source_asum_r.h"
#undef BASE_TYPE
}

float  gsl_blas_raw_scasum (size_t N, const gsl_complex_packed_array_float X, size_t incX)
{
#define BASE_TYPE float
#include "source_asum_c.h"
#undef BASE_TYPE
}

double gsl_blas_raw_dzasum (size_t N, const gsl_complex_packed_array X, size_t incX)
{
#define BASE_TYPE double
#include "source_asum_c.h"
#undef BASE_TYPE
}


CBLAS_INDEX gsl_blas_raw_isamax (size_t N, const float X[], size_t incX)
{
#define BASE_TYPE float
#include "source_iamax_r.h"
#undef BASE_TYPE
}

CBLAS_INDEX gsl_blas_raw_idamax (size_t N, const double X[], size_t incX)
{
#define BASE_TYPE double
#include "source_iamax_r.h"
#undef BASE_TYPE
}

CBLAS_INDEX gsl_blas_raw_icamax (size_t N, const gsl_complex_packed_array_float X, size_t incX)
{
#define BASE_TYPE float
#include "source_iamax_c.h"
#undef BASE_TYPE
}

CBLAS_INDEX gsl_blas_raw_izamax (size_t N, const gsl_complex_packed_array X, size_t incX)
{
#define BASE_TYPE double
#include "source_iamax_c.h"
#undef BASE_TYPE
}


void gsl_blas_raw_sswap (size_t N,
                         float X[], size_t incX,
                         float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_swap_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dswap (size_t N,
                         double X[], size_t incX,
                         double Y[], size_t incY)
{
#define BASE_TYPE double
#include "source_swap_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_cswap (size_t N,
                         gsl_complex_packed_array_float X, size_t incX,
                         gsl_complex_packed_array_float Y, size_t incY)
{
#define BASE_TYPE float
#include "source_swap_c.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zswap (size_t N,
                         gsl_complex_packed_array X, size_t incX,
                         gsl_complex_packed_array Y, size_t incY)
{
#define BASE_TYPE double
#include "source_swap_c.h"
#undef BASE_TYPE
}

void gsl_blas_raw_scopy (size_t N,
                         const float X[], size_t incX,
                         float Y[], size_t incY)
{
#include "source_copy_r.h"
}

void gsl_blas_raw_dcopy (size_t N,
                         const double X[], size_t incX,
                         double Y[], size_t incY)
{
#include "source_copy_r.h"
}

void gsl_blas_raw_ccopy (size_t N,
                         const gsl_complex_packed_array_float X, size_t incX,
                         gsl_complex_packed_array_float Y, size_t incY)
{
#include "source_copy_c.h"
}

void gsl_blas_raw_zcopy (size_t N,
                         const gsl_complex_packed_array X, size_t incX,
                         gsl_complex_packed_array Y, size_t incY)
{
#include "source_copy_c.h"
}

void gsl_blas_raw_saxpy (size_t N,
                         float alpha,
                         const float X[], size_t incX,
                         float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_axpy_r.h"
#undef BASE_TYPE
}


void gsl_blas_raw_daxpy (size_t N,
                         double alpha,
                         const double X[], size_t incX, 
                         double Y[], size_t incY)
{
#define BASE_TYPE double
#include "source_axpy_r.h"
#undef BASE_TYPE
}


void gsl_blas_raw_caxpy (size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float X, size_t incX,
                         gsl_complex_packed_array_float Y, size_t incY)
{
#define BASE_TYPE float
#include "source_axpy_c.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zaxpy (size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array X, size_t incX,
                         gsl_complex_packed_array Y, size_t incY)
{
#define BASE_TYPE double
#include "source_axpy_c.h"
#undef BASE_TYPE
}


void gsl_blas_raw_srotg (float * a, float * b, float * c, float * s)
{
#define BASE_TYPE float
#include "source_rotg.h"
#undef BASE_TYPE
}

void gsl_blas_raw_drotg (double a[], double b[], double c[], double s[])
{
#define BASE_TYPE double
#include "source_rotg.h"
#undef BASE_TYPE
}

void gsl_blas_raw_srotmg (float * d1, float * d2,
                          float * b1, float b2, float P[])
{
#define BASE_TYPE float
#include "source_rotmg.h"
#undef BASE_TYPE
}

void gsl_blas_raw_drotmg (double * d1, double * d2,
                          double * b1, double b2, double P[])
{
#define BASE_TYPE double
#include "source_rotmg.h"
#undef BASE_TYPE
}


void gsl_blas_raw_srot (size_t N,
                        float X[], size_t incX,
                        float Y[], size_t incY,
                        const float c, const float s)
{
#define BASE_TYPE float
#include "source_rot.h"
#undef BASE_TYPE
}
void gsl_blas_raw_drot (size_t N,
                        double X[], size_t incX,
                        double Y[], size_t incY,
                        const double c, const double s)
{
#define BASE_TYPE double
#include "source_rot.h"
#undef BASE_TYPE
}

void gsl_blas_raw_srotm (size_t N,
                         float X[], size_t incX,
                         float Y[], size_t incY,
                         const float P[])
{
#define BASE_TYPE float
#include "source_rotm.h"
#undef BASE_TYPE
}

void gsl_blas_raw_drotm (size_t N,
                         double X[], size_t incX,
                         double Y[], size_t incY,
                         const double P[])
{
#define BASE_TYPE double
#include "source_rotm.h"
#undef BASE_TYPE
}


void gsl_blas_raw_sscal  (size_t N, float  alpha, float  X[], size_t incX)
{
#include "source_scal_r.h"
}
void gsl_blas_raw_dscal  (size_t N, double alpha, double X[], size_t incX)
{
#include "source_scal_r.h"
}

void gsl_blas_raw_cscal  (size_t N, const gsl_complex_packed_float alpha, gsl_complex_packed_array_float X, size_t incX)
{
#define BASE_TYPE float
#include "source_scal_c.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zscal  (size_t N, const gsl_complex_packed alpha, gsl_complex_packed_array X, size_t incX)
{
#define BASE_TYPE double
#include "source_scal_c.h"
#undef BASE_TYPE
}

void gsl_blas_raw_csscal (size_t N, float  alpha, gsl_complex_packed_array_float X, size_t incX)
{
#include "source_scal_c_s.h"
}

void gsl_blas_raw_zdscal (size_t N, double alpha, gsl_complex_packed_array X, size_t incX)
{
#include "source_scal_c_s.h"
}


/*
 * ===========================================================================
 * level 2 BLAS
 * ===========================================================================
 */


/* GEMV */

void gsl_blas_raw_sgemv(CBLAS_TRANSPOSE TransA,
                        size_t M, size_t N,
                        float alpha,
                        const float A[], int lda,
                        const float X[], size_t incX,
                        float beta,
                        float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_gemv.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dgemv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
#define BASE_TYPE double
#include "source_gemv.h"
#undef BASE_TYPE
}

void gsl_blas_raw_cgemv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
#define BASE_TYPE float
#include "source_gemv_c.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zgemv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
#define BASE_TYPE double
#include "source_gemv_c.h"
#undef BASE_TYPE
}

/* GBMV */

void gsl_blas_raw_sgbmv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         float alpha,
                         const float A[], int lda,
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_gbmv_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dgbmv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
#define BASE_TYPE double
#include "source_gbmv_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_cgbmv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
#define BASE_TYPE float
#include "source_gbmv_c.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zgbmv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N, size_t KL, size_t KU,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
#define BASE_TYPE double
#include "source_gbmv_c.h"
#undef BASE_TYPE
}


/* TRMV */

void gsl_blas_raw_strmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA,
                         CBLAS_DIAG Diag,
                         size_t N,
                         const float A[], int lda,
                         float X[], size_t incX)
{
  gsl_blas_raw_stbmv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}

void gsl_blas_raw_dtrmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const double A[], int lda,
                         double X[], size_t incX)
{
  gsl_blas_raw_dtbmv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}

void gsl_blas_raw_ctrmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const gsl_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX)
{
  gsl_blas_raw_ctbmv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}

void gsl_blas_raw_ztrmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, size_t incX)
{
  gsl_blas_raw_ztbmv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}


/* TBMV */

void gsl_blas_raw_stbmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N, size_t K,
                         const float A[], int lda,
                         float X[], size_t incX)
{
#define BASE_TYPE float
#include "source_tbmv_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dtbmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N, size_t K,
                         const double A[], int lda,
                         double X[], size_t incX)
{
#define BASE_TYPE double
#include "source_tbmv_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_ctbmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N, size_t K,
                         const gsl_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX)
{
/* FIXME: copy above when debugged */
}

void gsl_blas_raw_ztbmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N, size_t K,
                         const gsl_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX)
{
/* FIXME: copy above when debugged */
}


/* TPMV */

void gsl_blas_raw_stpmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const float Ap[],
                         float X[], size_t incX)
{
#define BASE_TYPE float
#include "source_tpmv_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dtpmv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const double Ap[],
                         double X[], size_t incX)
{
#define BASE_TYPE double
#include "source_tpmv_r.h"
#undef BASE_TYPE
}

void gsl_blas_raw_ctpmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float Ap,
                               gsl_complex_packed_array_float X, size_t incX)
{
/* FIXME: copy above when debugged */
}

void gsl_blas_raw_ztpmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array Ap,
                               gsl_complex_packed_array X, size_t incX)
{
/* FIXME: copy above when debugged */
}

/* TRSV */

void gsl_blas_raw_strsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const float A[], int lda,
                         float X[], size_t incX)
{
  gsl_blas_raw_stbsv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}

void gsl_blas_raw_dtrsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const double A[], int lda,
                         double X[], size_t incX)
{
  gsl_blas_raw_dtbsv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}

void gsl_blas_raw_ctrsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const gsl_complex_packed_array_float A, int lda,
                         gsl_complex_packed_array_float X, size_t incX)
{
  gsl_blas_raw_ctbsv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}

void gsl_blas_raw_ztrsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const gsl_complex_packed_array A, int lda,
                         gsl_complex_packed_array X, size_t incX)
{
  gsl_blas_raw_ztbsv (Uplo,
                      TransA, Diag,
                      N, N-1,
                      A, lda,
                      X, incX);
}


/* TBSV */

void gsl_blas_raw_stbsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N, size_t K,
                         const float A[], int lda,
                         float X[], size_t incX)
{
#define BASE_TYPE float
#define MATRIX_VAR_NAME A
#include "source_tbsv_r.h"
#undef MATRIX_VAR_NAME
#undef BASE_TYPE
}

void gsl_blas_raw_dtbsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N, size_t K,
                         const double A[], int lda,
                         double X[], size_t incX)
{
#define BASE_TYPE double
#define MATRIX_VAR_NAME A
#include "source_tbsv_r.h"
#undef MATRIX_VAR_NAME
#undef BASE_TYPE
}

void gsl_blas_raw_ctbsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, size_t incX)
{
/* FIXME: copy above */
}

void gsl_blas_raw_ztbsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, size_t incX)
{
/* FIXME: copy above */
}


/* TPSV */

void gsl_blas_raw_stpsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const float Ap[],
                         float X[], size_t incX)
{
#define BASE_TYPE float
#define MATRIX_VAR_NAME Ap
#include "source_tpsv_r.h"
#undef MATRIX_VAR_NAME
#undef BASE_TYPE
}

void gsl_blas_raw_dtpsv (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                         size_t N,
                         const double Ap[],
                         double X[], size_t incX)
{
#define BASE_TYPE double
#define MATRIX_VAR_NAME Ap
#include "source_tpsv_r.h"
#undef MATRIX_VAR_NAME
#undef BASE_TYPE
}

void gsl_blas_raw_ctpsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float Ap,
                               gsl_complex_packed_array_float X, size_t incX)
{
/* FIXME: copy above */
}

void gsl_blas_raw_ztpsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array Ap,
                               gsl_complex_packed_array X, size_t incX)
{
/* FIXME: copy above */
}


/* SYMV */

void gsl_blas_raw_ssymv (CBLAS_UPLO Uplo,
                         size_t N,
                         float alpha,
                         const float A[], int lda,
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_symv.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dsymv (CBLAS_UPLO Uplo,
                               size_t N,
                               double alpha,
                               const double A[], int lda,
                               const double X[], size_t incX,
                               double beta,
                               double Y[], size_t incY)
{
#define BASE_TYPE double
#include "source_symv.h"
#undef BASE_TYPE
}


/* SBMV */

void gsl_blas_raw_ssbmv (CBLAS_UPLO Uplo,
                         size_t N, size_t K,
                         float alpha,
                         const float A[], int lda, 
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_sbmv.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dsbmv (CBLAS_UPLO Uplo,
                         size_t N, size_t K,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
#define BASE_TYPE double
#include "source_sbmv.h"
#undef BASE_TYPE
}

/* SPMV */

void gsl_blas_raw_sspmv (CBLAS_UPLO Uplo,
                         size_t N,
                         float alpha,
                         const float Ap[],
                         const float X[], size_t incX,
                         float beta,
                         float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_spmv.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dspmv (CBLAS_UPLO Uplo,
                         size_t N,
                         double alpha,
                         const double Ap[],
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
#define BASE_TYPE double
#include "source_spmv.h"
#undef BASE_TYPE
}


/* GER */

void gsl_blas_raw_sger (size_t M, size_t N,
                        float alpha,
                        const float X[], size_t incX,
                        const float Y[], size_t incY,
                        float A[], int lda)
{
#define BASE_TYPE float
#include "source_ger.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dger (size_t M, size_t N,
                        double alpha,
                        const double X[], size_t incX,
                        const double Y[], size_t incY,
                        double A[], int lda)
{
#define BASE_TYPE double
#include "source_ger.h"
#undef BASE_TYPE
}


/* SYR */

void gsl_blas_raw_ssyr (CBLAS_UPLO Uplo,
                        size_t N,
                        float alpha,
                        const float X[], size_t incX,
                        float A[], int lda)
{
#define BASE_TYPE float
#include "source_syr.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dsyr (CBLAS_UPLO Uplo,
                        size_t N,
                        double alpha,
                        const double X[], size_t incX,
                        double A[], int lda)
{
#define BASE_TYPE double
#include "source_syr.h"
#undef BASE_TYPE
}


/* SPR */

void gsl_blas_raw_sspr (CBLAS_UPLO Uplo,
                        size_t N,
                        float alpha,
                        const float X[], size_t incX,
                        float Ap[])
{
#define BASE_TYPE float
#include "source_spr.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dspr (CBLAS_UPLO Uplo,
                        size_t N,
                        double alpha,
                        const double X[], size_t incX,
                        double Ap[])
{
#define BASE_TYPE double
#include "source_spr.h"
#undef BASE_TYPE
}


/* SYR2 */

void gsl_blas_raw_ssyr2 (CBLAS_UPLO Uplo,
                         size_t N,
                         float alpha,
                         const float X[], size_t incX, 
                         const float Y[], size_t incY,
                         float A[], int lda)
{
#define BASE_TYPE float
#include "source_syr2.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dsyr2 (CBLAS_UPLO Uplo,
                               size_t N,
                               double alpha,
                               const double X[], size_t incX,
                               const double Y[], size_t incY,
                               double A[], int lda)
{
#define BASE_TYPE double
#include "source_syr2.h"
#undef BASE_TYPE
}


/* SPR2 */

void gsl_blas_raw_sspr2 (CBLAS_UPLO Uplo,
                         size_t N,
                         float alpha,
                         const float X[], size_t incX,
                         const float Y[], size_t incY,
                         float A[])
{
#define BASE_TYPE double
#include "source_spr2.h"
#undef BASE_TYPE
}

void gsl_blas_raw_dspr2 (CBLAS_UPLO Uplo,
                         size_t N,
                         double alpha,
                         const double X[], size_t incX,
                         const double Y[], size_t incY,
                         double A[])
{
#define BASE_TYPE double
#include "source_spr2.h"
#undef BASE_TYPE
}


/* HEMV */

void gsl_blas_raw_chemv (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
  gsl_blas_raw_chbmv (Uplo,
                      N, N-1,
                      alpha,
                      A, lda,
                      X, incX,
                      beta,
                      Y, incY);
}

void gsl_blas_raw_zhemv (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
  gsl_blas_raw_zhbmv (Uplo,
                      N, N-1,
                      alpha,
                      A, lda,
                      X, incX,
                      beta,
                      Y, incY);
}


/* HBMV */

void gsl_blas_raw_chbmv (CBLAS_UPLO Uplo,
                         size_t N, size_t K,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
#define BASE_TYPE float
#include "source_hbmv.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zhbmv (CBLAS_UPLO Uplo,
                         size_t N, size_t K,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
#define BASE_TYPE double
#include "source_hbmv.h"
#undef BASE_TYPE
}


/* HPMV */

void gsl_blas_raw_chpmv (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float Ap,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float Y, size_t incY)
{
#define BASE_TYPE float
#include "source_hpmv.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zhpmv (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array Ap,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
#define BASE_TYPE double
#include "source_hpmv.h"
#undef BASE_TYPE
}


/* GERU */

void gsl_blas_raw_cgeru (size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float A, int lda)
{
#define BASE_TYPE float
#include "source_geru.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zgeru (size_t M, size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array A, int lda)
{
#define BASE_TYPE double
#include "source_geru.h"
#undef BASE_TYPE
}


/* GERC */

void gsl_blas_raw_cgerc (size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float A, int lda)
{
#define BASE_TYPE float
#include "source_gerc.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zgerc (size_t M, size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array A, int lda)
{
#define BASE_TYPE double
#include "source_gerc.h"
#undef BASE_TYPE
}

/* HER */

void gsl_blas_raw_cher (CBLAS_UPLO Uplo,
                        size_t N,
                        float alpha,
                        const gsl_complex_packed_array_float X, size_t incX,
                        gsl_complex_packed_array_float A, int lda)
{
#define BASE_TYPE float
#include "source_her.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zher (CBLAS_UPLO Uplo,
                        size_t N,
                        double alpha,
                        const gsl_complex_packed_array X, size_t incX,
                        gsl_complex_packed_array A, int lda)
{
#define BASE_TYPE double
#include "source_her.h"
#undef BASE_TYPE
}


/* HPR */

void gsl_blas_raw_chpr (CBLAS_UPLO Uplo,
                        size_t N,
                        float alpha,
                        const gsl_complex_packed_array_float X, size_t incX,
                        gsl_complex_packed_array_float Ap)
{
#define BASE_TYPE float
#include "source_hpr.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zhpr (CBLAS_UPLO Uplo,
                        size_t N,
                        double alpha,
                        const gsl_complex_packed_array X, size_t incX,
                        gsl_complex_packed_array Ap)
{
#define BASE_TYPE double
#include "source_hpr.h"
#undef BASE_TYPE
}


/* HER2 */

void gsl_blas_raw_cher2 (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float A, int lda)
{
#define BASE_TYPE float
#include "source_her2.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zher2 (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array A, int lda)
{
#define BASE_TYPE double
#include "source_her2.h"
#undef BASE_TYPE
}


/* HPR2 */

void gsl_blas_raw_chpr2 (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float X, size_t incX,
                         const gsl_complex_packed_array_float Y, size_t incY,
                         gsl_complex_packed_array_float Ap)
{
#define BASE_TYPE float
#include "source_hpr2.h"
#undef BASE_TYPE
}

void gsl_blas_raw_zhpr2 (CBLAS_UPLO Uplo,
                         size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed_array Y, size_t incY,
                         gsl_complex_packed_array Ap)
{
#define BASE_TYPE double
#include "source_hpr2.h"
#undef BASE_TYPE
}


/*
 * ===========================================================================
 * level 3 BLAS
 * ===========================================================================
 */

/* GEMM */

void gsl_blas_raw_sgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
                         size_t M, size_t N, size_t K,
                         float alpha,
                         const float A[], int lda,
                         const float B[], int ldb,
                         float beta,
                         float C[], int ldc)
{
}

void gsl_blas_raw_dgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
                         size_t M, size_t N, size_t K,
                         double alpha,
                         const double A[], int lda,
                         const double B[], int ldb,
                         double beta,
                         double C[], int ldc)
{
}

void gsl_blas_raw_cgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
                         size_t M, size_t N, size_t K,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float A, int lda, 
                         const gsl_complex_packed_array_float B, int ldb,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
}

void gsl_blas_raw_zgemm (CBLAS_TRANSPOSE TransA,
                         CBLAS_TRANSPOSE TransB,
                         size_t M, size_t N, size_t K,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array B, int ldb,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array C, int ldc)
{
}


/* SYMM */

void gsl_blas_raw_ssymm (CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                         size_t M, size_t N,
                         float alpha,
                         const float A[], int lda,
                         const float B[], int ldb,
                         float beta,
                         float C[], int ldc)
{
}

void gsl_blas_raw_dsymm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
                         size_t M, size_t N,
                         double alpha,
                         const double A[], int lda,
                         const double B[], int ldb,
                         double beta,
                         double C[], int ldc)
{
}

void gsl_blas_raw_csymm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
                         size_t M, size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_array_float B, int ldb,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
}

void gsl_blas_raw_zsymm (CBLAS_SIDE Side,
                         CBLAS_UPLO Uplo,
                         size_t M, size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array B, int ldb,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array C, int ldc)
{
}


/* SYRK */

void gsl_blas_raw_ssyrk (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                         size_t N, size_t K,
                         float alpha,
                         const float A[], int lda,
                         float beta,
                         float C[], int ldc)
{
}

void gsl_blas_raw_dsyrk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
                         size_t N, size_t K,
                         double alpha,
                         const double A[], int lda,
                         double beta,
                         double C[], int ldc)
{
}

void gsl_blas_raw_csyrk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
                         size_t N, size_t K,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float A, int lda,
                         const gsl_complex_packed_float beta,
                         gsl_complex_packed_array_float C, int ldc)
{
}

void gsl_blas_raw_zsyrk (CBLAS_UPLO Uplo,
                         CBLAS_TRANSPOSE Trans,
                         size_t N, size_t K,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array C, int ldc)
{
}


/* SYR2K */

void gsl_blas_raw_ssyr2k (CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
                          size_t N, size_t K,
                          float alpha,
                          const float A[], int lda,
                          const float B[], int ldb,
                          float beta,
                          float C[], int ldc)
{
}

void gsl_blas_raw_dsyr2k (CBLAS_UPLO Uplo,
                          CBLAS_TRANSPOSE Trans,
                          size_t N, size_t K,
                          double alpha,
                          const double A[], int lda,
                          const double B[], int ldb,
                          double beta,
                          double C[], int ldc)
{
}

void gsl_blas_raw_csyr2k (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                const gsl_complex_packed_float alpha,
                                const gsl_complex_packed_array_float A, int lda,
                                const gsl_complex_packed_array_float B, int ldb,
                                const gsl_complex_packed_float beta,
                                gsl_complex_packed_array_float C, int ldc)
{
}

void gsl_blas_raw_zsyr2k (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                const gsl_complex_packed alpha,
                                const gsl_complex_packed_array A, int lda,
                                const gsl_complex_packed_array B, int ldb,
                                const gsl_complex_packed beta,
                                gsl_complex_packed_array C, int ldc)
{
}


/* TRMM */

void gsl_blas_raw_strmm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               float alpha,
                               const float A[], int lda,
                               float B[], int ldb)
{
}

void gsl_blas_raw_dtrmm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               double alpha,
                               const double A[], int lda,
                               double B[], int ldb)
{
}

void gsl_blas_raw_ctrmm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float B, int ldb)
{
}

void gsl_blas_raw_ztrmm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array B, int ldb)
{
}


/* TRSM */

void gsl_blas_raw_strsm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               float alpha,
                               const float A[], int lda,
                               float B[], int ldb)
{
}

void gsl_blas_raw_dtrsm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               double alpha,
                               const double A[], int lda,
                               double B[], int ldb)
{
}

void gsl_blas_raw_ctrsm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float B, int ldb)
{
}

void gsl_blas_raw_ztrsm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array B, int ldb)
{
}


/* HEMM */

void gsl_blas_raw_chemm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float B, int ldb,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float C, int ldc)
{
}

void gsl_blas_raw_zhemm (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo,
                               size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array B, int ldb,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array C, int ldc)
{
}


/* HERK */

void gsl_blas_raw_cherk (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               float beta,
                               gsl_complex_packed_array_float C, int ldc)
{
}

void gsl_blas_raw_zherk (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               double alpha,
                               const gsl_complex_packed_array A, int lda,
                               double beta,
                               gsl_complex_packed_array C, int ldc)
{
}


/* HER2K */

void gsl_blas_raw_cher2k (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                const gsl_complex_packed_float alpha,
                                const gsl_complex_packed_array_float A, int lda,
                                const gsl_complex_packed_array_float B, int ldb,
                                float beta,
                                gsl_complex_packed_array_float C, int ldc)
{
}


void gsl_blas_raw_zher2k (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                const gsl_complex_packed alpha,
                                const gsl_complex_packed_array A, int lda,
                                const gsl_complex_packed_array B, int ldb,
                                double beta,
                                gsl_complex_packed_array C, int ldc)
{
}
