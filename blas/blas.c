/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include "gsl_blas_types.h"
#include "gsl_blas_raw.h"
#include "gsl_blas.h"


/* ========================================================================
 * Level 1
 * ========================================================================
 */

int gsl_blas_sdsdot (float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
		     float * result
		     )
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  *result = gsl_blas_raw_sdsdot(N, alpha, X->data, X->stride, Y->data, Y->stride);
  return status;
}

int gsl_blas_dsdot (const gsl_vector_float * X,
                    const gsl_vector_float * Y,
		    double * result
		    )
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  *result = gsl_blas_raw_dsdot(N, X->data, X->stride, Y->data, Y->stride);
  return status;
}

int gsl_blas_sdot (const gsl_vector_float * X,
                   const gsl_vector_float * Y,
		   float * result
		   )
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  *result = gsl_blas_raw_sdot(N, X->data, X->stride, Y->data, Y->stride);
  return status;
}

int gsl_blas_ddot (const gsl_vector * X,
                   const gsl_vector * Y,
		   double * result
		   )
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  *result = gsl_blas_raw_ddot(N, X->data, X->stride, Y->data, Y->stride);
  return status;
}


/*
 * Functions having prefixes Z and C only
 */
int  gsl_blas_cdotu (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotu)
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_cdotu(N, X->data, X->stride, Y->data, Y->stride, dotu);
  return status;
}


int  gsl_blas_cdotc (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotc)
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_cdotc(N, X->data, X->stride, Y->data, Y->stride, dotc);
  return status;
}


int  gsl_blas_zdotu (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotu)
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_zdotu(N, X->data, X->stride, Y->data, Y->stride, dotu);
  return status;
}


int  gsl_blas_zdotc (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotc)
{
  CBLAS_INDEX N = MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_zdotc(N, X->data, X->stride, Y->data, Y->stride, dotc);
  return status;
}


/*
 * Functions having prefixes S D SC DZ
 */
float  gsl_blas_snrm2 (const gsl_vector_float * X)
{
  return gsl_blas_raw_snrm2(X->size, X->data, X->stride);
}

float  gsl_blas_sasum (const gsl_vector_float * X)
{
  return gsl_blas_raw_sasum(X->size, X->data, X->stride);
}

double gsl_blas_dnrm2 (const gsl_vector * X)
{
  return gsl_blas_raw_dnrm2(X->size, X->data, X->stride);
}

double gsl_blas_dasum (const gsl_vector * X)
{
  return gsl_blas_raw_dasum(X->size, X->data, X->stride);
}

float  gsl_blas_scnrm2 (const gsl_vector_complex * X)
{
  return gsl_blas_raw_scnrm2(X->size, X->data, X->stride);
}

float  gsl_blas_scasum (const gsl_vector_complex * X)
{
  return gsl_blas_raw_scasum(X->size, X->data, X->stride);
}

double gsl_blas_dznrm2 (const gsl_vector_complex * X)
{
  return gsl_blas_raw_dznrm2(X->size, X->data, X->stride);
}

double gsl_blas_dzasum (const gsl_vector_complex * X)
{
  return gsl_blas_raw_dzasum(X->size, X->data, X->stride);
}


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX gsl_blas_isamax (const gsl_vector_float * X)
{
  return gsl_blas_raw_isamax(X->size, X->data, X->stride);
}

CBLAS_INDEX gsl_blas_idamax (const gsl_vector * X)
{
  return gsl_blas_raw_idamax(X->size, X->data, X->stride);
}

CBLAS_INDEX gsl_blas_icamax (const gsl_vector_complex * X)
{
  return gsl_blas_raw_icamax(X->size, X->data, X->stride);
}

CBLAS_INDEX gsl_blas_izamax (const gsl_vector_complex * X)
{
  return gsl_blas_raw_izamax(X->size, X->data, X->stride);
}



/*
 * Routines with standard 4 prefixes (s, d, c, z)
 */
int  gsl_blas_sswap (gsl_vector_float * X,
                     gsl_vector_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_sswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_scopy (const gsl_vector_float * X,
                     gsl_vector_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_scopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_saxpy (float alpha, const gsl_vector_float * X, gsl_vector_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_saxpy(X->size, alpha, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_dswap (gsl_vector * X, gsl_vector * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_dswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_dcopy (const gsl_vector * X, gsl_vector * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_dcopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_daxpy (double alpha, const gsl_vector * X, gsl_vector * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_daxpy(X->size, alpha, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_cswap (gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_cswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_ccopy (const gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_ccopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_caxpy (const gsl_complex * alpha, const gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_caxpy(X->size, alpha, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_zswap (gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_zswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_zcopy (const gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_zcopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_zaxpy (const gsl_complex * alpha, const gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_zaxpy(X->size, alpha, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/*
 * Routines with S and D prefix only
 */
int
gsl_blas_srotg (float a[], float b[], float c[], float s[])
{
}


int
gsl_blas_srotmg (float d1[], float d2[], float b1[], float b2, float P[])
{
}


int
gsl_blas_srot (gsl_vector_float * X, gsl_vector_float * Y, float c, float s)
{
}


int
gsl_blas_srotm (gsl_vector_float * X, gsl_vector_float * Y, const float P[])
{
}


int
gsl_blas_drotg (double a[], double b[], double c[], double s[])
{
}


int
gsl_blas_drotmg (double d1[], double d2[], double b1[], double b2, double P[])
{
}


int
gsl_blas_drot (gsl_vector * X, gsl_vector * Y, const double c, const double s)
{
}


int
gsl_blas_drotm (gsl_vector * X, gsl_vector * Y, const double P[])
{
}


/*
 * Routines with S D C Z CS and ZD prefixes
 */
void gsl_blas_sscal (float  alpha, gsl_vector_float  * X)
{
  gsl_blas_raw_sscal(X->size, alpha, X->data, X->stride);
}

void gsl_blas_dscal (double alpha, gsl_vector * X)
{
  gsl_blas_raw_dscal(X->size, alpha, X->data, X->stride);
}

void gsl_blas_cscal  (const gsl_complex * alpha, gsl_vector_complex * X)
{
  gsl_blas_raw_cscal(X->size, alpha, X->data, X->stride);
}

void gsl_blas_zscal  (const gsl_complex * alpha, gsl_vector_complex * X)
{
  gsl_blas_raw_zscal(X->size, alpha, X->data, X->stride);
}

void gsl_blas_csscal (float  alpha, gsl_vector_complex * X)
{
  gsl_blas_raw_csscal(X->size, alpha, X->data, X->stride);
}

void gsl_blas_zdscal (double alpha, gsl_vector_complex * X)
{
  gsl_blas_raw_zdscal(X->size, alpha, X->data, X->stride);
}




/* ===========================================================================
 * Level 2
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
int  gsl_blas_sgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
}


int  gsl_blas_sgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
		     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
}


int  gsl_blas_strmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
}


int  gsl_blas_stbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
}


int  gsl_blas_stpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const float Ap[],
                     gsl_vector_float * X)
{
}


int  gsl_blas_strsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
}


int  gsl_blas_stbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
}


int  gsl_blas_stpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const float Ap[],
                     gsl_vector_float * X)
{
}


int  gsl_blas_dgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
}


int  gsl_blas_dgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
}


int  gsl_blas_dtrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
}


int  gsl_blas_dtbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
}


int  gsl_blas_dtpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const double Ap[],
                     gsl_vector * X)
{
}


int  gsl_blas_dtrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
}


int  gsl_blas_dtbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
}


int  gsl_blas_dtpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const double Ap[],
                     gsl_vector * X)
{
}


int  gsl_blas_cgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_cgbmv (CBLAS_ORDER order,
                     CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_ctrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix * A,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ctbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ctpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void *Ap,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ctrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ctbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ctpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void * Ap,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_zgemv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_zgbmv (CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_ztrmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ztbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ztpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void *Ap,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ztrsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex *X)
{
}


int  gsl_blas_ztbsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
}


int  gsl_blas_ztpsv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void *Ap,
                     gsl_vector_complex * X)
{
}


/*
 * Routines with S and D prefixes only
 */
int  gsl_blas_ssymv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
}


int  gsl_blas_ssbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
}


int  gsl_blas_sspmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     float alpha,
                     const float Ap[],
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
}


int  gsl_blas_sger (CBLAS_ORDER order,
                    float alpha,
                    const gsl_vector_float * X,
                    const gsl_vector_float * Y,
                    gsl_matrix_float * A)
{
}


int  gsl_blas_ssyr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    gsl_matrix_float * A)
{
}


int  gsl_blas_sspr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    float Ap[])
{
}


int  gsl_blas_ssyr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     gsl_matrix_float * A)
{
}


int  gsl_blas_sspr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     gsl_vector_float * A)
{
}


int  gsl_blas_dsymv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
}


int  gsl_blas_dsbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
}


int  gsl_blas_dspmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     double alpha,
                     const double Ap[],
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
}


int  gsl_blas_dger (CBLAS_ORDER order,
                    double alpha,
                    const gsl_vector * X,
                    const gsl_vector * Y,
                    gsl_matrix * A)
{
}


int  gsl_blas_dsyr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector * X,
                    gsl_matrix * A)
{
}


int  gsl_blas_dspr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector * X,
                    double Ap[])
{
}


int  gsl_blas_dsyr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     gsl_matrix * A)
{
}


int  gsl_blas_dspr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     double Ap[])
{
}


/*
 * Routines with C and Z prefixes only
 */
int  gsl_blas_chemv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_chbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_chpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
		     const void * Ap,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_cgeru (CBLAS_ORDER order,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


int  gsl_blas_cgerc (CBLAS_ORDER order,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


int  gsl_blas_cher (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_complex * X,
                    gsl_matrix_complex * A)
{
}


int  gsl_blas_chpr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_complex * X,
                    void * Ap)
{
}


int  gsl_blas_cher2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


int  gsl_blas_chpr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     void * Ap)
{
}


int  gsl_blas_zhemv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_zhbmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_zhpmv (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const void * Ap,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


int  gsl_blas_zgeru (CBLAS_ORDER order,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


int  gsl_blas_zgerc (CBLAS_ORDER order,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


int  gsl_blas_zher (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    gsl_matrix_complex * A)
{
}


int  gsl_blas_zhpr (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    void * Ap)
{
}


int  gsl_blas_zher2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


int  gsl_blas_zhpr2 (CBLAS_ORDER order, CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     void *Ap)
{
}


/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/*
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
int  gsl_blas_sgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C)
{
}


int  gsl_blas_ssymm (CBLAS_ORDER Order, CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C)
{
}


int  gsl_blas_ssyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     float beta,
                     gsl_matrix_float * C)
{
}


int  gsl_blas_ssyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                      int K,
                      float alpha,
                      const gsl_matrix_float * A,
                      const gsl_matrix_float * B,
                      float beta,
                      gsl_matrix_float * C)
{
}


int  gsl_blas_strmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B)
{
}


int  gsl_blas_strsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B)
{
}


int  gsl_blas_dgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C)
{
}


int  gsl_blas_dsymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C)
{
}


int  gsl_blas_dsyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     double beta,
                     gsl_matrix * C)
{
}


int  gsl_blas_dsyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      double alpha,
                      const  gsl_matrix * A,
                      const  gsl_matrix * B,
                      double beta,
                      gsl_matrix * C)
{
}


int  gsl_blas_dtrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B)
{
}


int  gsl_blas_dtrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B)
{
}


int  gsl_blas_cgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_csymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_csyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_csyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      const gsl_complex * beta,
                      gsl_matrix_complex * C)
{
}


int  gsl_blas_ctrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


int  gsl_blas_ctrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


int  gsl_blas_zgemm (CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_zsymm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_zsyrk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_zsyr2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      const gsl_complex * beta,
                      gsl_matrix_complex *C)
{
}


int  gsl_blas_ztrmm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


int  gsl_blas_ztrsm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


/*
 * Routines with prefixes C and Z only
 */
int  gsl_blas_chemm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_cherk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_complex * A,
                     float beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_cher2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      float beta,
                      gsl_matrix_complex * C)
{
}


int  gsl_blas_zhemm (CBLAS_ORDER Order, CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_zherk (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     double alpha,
                     const gsl_matrix_complex * A,
                     double beta,
                     gsl_matrix_complex * C)
{
}


int  gsl_blas_zher2k (CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      double beta,
                      gsl_matrix_complex * C)
{
}

