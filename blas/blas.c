/*
 * Author:  G. Jungman
 * RCS:     $Id$
 */
/* GSL implementation of BLAS operations.
 * Conforms to gsl_blas interface.
 * Note that GSL native storage is row-major, so
 * we implement in terms of gsl_blas_raw.
 */
#include <gsl_math.h>
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
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
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
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
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
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
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
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  *result = gsl_blas_raw_ddot(N, X->data, X->stride, Y->data, Y->stride);
  return status;
}


int  gsl_blas_cdotu (const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_complex_float * dotu)
{
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_cdotu(N, X->data, X->stride, Y->data, Y->stride, dotu->dat);
  return status;
}


int  gsl_blas_cdotc (const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_complex_float * dotc)
{
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_cdotc(N, X->data, X->stride, Y->data, Y->stride, dotc->dat);
  return status;
}


int  gsl_blas_zdotu (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotu)
{
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_zdotu(N, X->data, X->stride, Y->data, Y->stride, dotu->dat);
  return status;
}


int  gsl_blas_zdotc (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotc)
{
  CBLAS_INDEX N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EINVAL;
  gsl_blas_raw_zdotc(N, X->data, X->stride, Y->data, Y->stride, dotc->dat);
  return status;
}


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

float  gsl_blas_scnrm2 (const gsl_vector_complex_float * X)
{
  return gsl_blas_raw_scnrm2(X->size, X->data, X->stride);
}

float  gsl_blas_scasum (const gsl_vector_complex_float * X)
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


CBLAS_INDEX gsl_blas_isamax (const gsl_vector_float * X)
{
  return gsl_blas_raw_isamax(X->size, X->data, X->stride);
}

CBLAS_INDEX gsl_blas_idamax (const gsl_vector * X)
{
  return gsl_blas_raw_idamax(X->size, X->data, X->stride);
}

CBLAS_INDEX gsl_blas_icamax (const gsl_vector_complex_float * X)
{
  return gsl_blas_raw_icamax(X->size, X->data, X->stride);
}

CBLAS_INDEX gsl_blas_izamax (const gsl_vector_complex * X)
{
  return gsl_blas_raw_izamax(X->size, X->data, X->stride);
}


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
gsl_blas_cswap (gsl_vector_complex_float * X, gsl_vector_complex_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_cswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_ccopy (const gsl_vector_complex_float * X, gsl_vector_complex_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_ccopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int
gsl_blas_caxpy (const gsl_complex_float * alpha,
                const gsl_vector_complex_float * X,
		gsl_vector_complex_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_caxpy(X->size, (float *)alpha->dat, X->data, X->stride, Y->data, Y->stride);
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
    gsl_blas_raw_zaxpy(X->size, (double *)alpha->dat, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


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


void gsl_blas_sscal (float  alpha, gsl_vector_float  * X)
{
  gsl_blas_raw_sscal(X->size, alpha, X->data, X->stride);
}

void gsl_blas_dscal (double alpha, gsl_vector * X)
{
  gsl_blas_raw_dscal(X->size, alpha, X->data, X->stride);
}

void gsl_blas_cscal  (const gsl_complex_float * alpha,
                      gsl_vector_complex_float * X)
{
  gsl_blas_raw_cscal(X->size, (float *)alpha->dat, X->data, X->stride);
}

void gsl_blas_zscal  (const gsl_complex * alpha, gsl_vector_complex * X)
{
  gsl_blas_raw_zscal(X->size, (double *)alpha->dat, X->data, X->stride);
}

void gsl_blas_csscal (float  alpha, gsl_vector_complex_float * X)
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

/* GEMV */

int  gsl_blas_sgemv (CBLAS_TRANSPOSE TransA,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;
  
  if(col_dim == X->size) {
    gsl_blas_raw_sgemv(TransA,
                       row_dim, col_dim,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_dgemv (CBLAS_TRANSPOSE TransA,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(col_dim == X->size) {
    gsl_blas_raw_dgemv(TransA,
                       row_dim, col_dim,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_cgemv (CBLAS_TRANSPOSE TransA,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;
  
  if(col_dim == X->size) {
    gsl_blas_raw_cgemv(TransA,
                       row_dim, col_dim,
		       (float *)alpha->dat,
		       A->data, tda,
		       X->data, X->stride,
		       (float *)beta->dat,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_zgemv (CBLAS_TRANSPOSE TransA,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;
  
  if(col_dim == X->size) {
    gsl_blas_raw_zgemv(TransA,
                       row_dim, col_dim,
		       (double *)alpha->dat,
		       A->data, tda,
		       X->data, X->stride,
		       (double *)beta->dat,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/* GBMV */

int  gsl_blas_sgbmv (CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
		     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;
  
  if(col_dim == X->size) {
    gsl_blas_raw_sgbmv(TransA,
                       row_dim, col_dim, KL, KU,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}

int  gsl_blas_dgbmv (CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(col_dim == X->size) {
    gsl_blas_raw_dgbmv(TransA,
                       row_dim, col_dim, KL, KU,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_cgbmv (CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(col_dim == X->size) {
    gsl_blas_raw_cgbmv(TransA,
                       row_dim, col_dim, KL, KU,
		       (float *)alpha->dat,
		       A->data, tda,
		       X->data, X->stride,
		       (float *)beta->dat,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_zgbmv (CBLAS_TRANSPOSE TransA,
                     int KL, int KU,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;
  
  if(col_dim == X->size) {
    gsl_blas_raw_zgbmv(TransA,
                       row_dim, col_dim, KL, KU,
		       (double *)alpha->dat,
		       A->data, tda,
		       X->data, X->stride,
		       (double *)beta->dat,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/* TRMV */

int  gsl_blas_strmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_strmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_dtrmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_dtrmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_ctrmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ctrmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_ztrmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ztrmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/* TBMV */

int  gsl_blas_stbmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_stbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_dtbmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_dtbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_ctbmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ctbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_ztbmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ztbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/* TPMV */

int  gsl_blas_stpmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const float Ap[],
                     gsl_vector_float * X)
{
}


int  gsl_blas_dtpmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const double Ap[],
                     gsl_vector * X)
{
}


int  gsl_blas_ctpmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void *Ap,
                     gsl_vector_complex_float * X)
{
}


int  gsl_blas_ztpmv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void *Ap,
                     gsl_vector_complex * X)
{
}


/* TRSV */

int  gsl_blas_strsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_strsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL; 
}


int  gsl_blas_dtrsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_dtrsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL; 
}


int  gsl_blas_ctrsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ctrsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL; 
}


int  gsl_blas_ztrsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex *X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ztrsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL; 
}


/* TBSV */

int  gsl_blas_stbsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_stbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_dtbsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_dtbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_ctbsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ctbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_ztbsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ztbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/* TPSV */

int  gsl_blas_stpsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const float Ap[],
                     gsl_vector_float * X)
{
}


int  gsl_blas_dtpsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const double Ap[],
                     gsl_vector * X)
{
}


int  gsl_blas_ctpsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void * Ap,
                     gsl_vector_complex_float * X)
{
}


int  gsl_blas_ztpsv (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                     const void *Ap,
                     gsl_vector_complex * X)
{
}


/* SYMV */

int  gsl_blas_ssymv (CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size1;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ssymv(Uplo,
                       col_dim,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_dsymv (CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size1;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_dsymv(Uplo,
                       col_dim,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/* SBMV */

int  gsl_blas_ssbmv (CBLAS_UPLO Uplo,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size1;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_ssbmv(Uplo,
                       col_dim,
		       K,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


int  gsl_blas_dsbmv (CBLAS_UPLO Uplo,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size1;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;

  if(col_dim == X->size) {
    gsl_blas_raw_dsbmv(Uplo,
                       col_dim,
		       K,
		       alpha,
		       A->data, tda,
		       X->data, X->stride,
		       beta,
		       Y->data, Y->stride
		       );
    return GSL_SUCCESS;
  }
  else
    return GSL_EINVAL;
}


/* SPMV */

int  gsl_blas_sspmv (CBLAS_UPLO Uplo,
                     float alpha,
                     const float Ap[],
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
}


int  gsl_blas_dspmv (CBLAS_UPLO Uplo,
                     double alpha,
                     const double Ap[],
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
}


/* GER */

int  gsl_blas_sger (float alpha,
                    const gsl_vector_float * X,
                    const gsl_vector_float * Y,
                    gsl_matrix_float * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->dim2 ;

  /* FIXME: check this */
  if(X->size != row_dim) return GSL_EINVAL;
  if(Y->size != col_dim) return GSL_EINVAL;

  gsl_blas_raw_sger(row_dim, col_dim,
		    alpha,
        	    X->data, X->stride,
        	    Y->data, Y->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


int  gsl_blas_dger (double alpha,
                    const gsl_vector * X,
                    const gsl_vector * Y,
                    gsl_matrix * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->dim2 ;

  /* FIXME: check this */
  if(X->size != row_dim) return GSL_EINVAL;
  if(Y->size != col_dim) return GSL_EINVAL;

  gsl_blas_raw_dger(row_dim, col_dim,
		    alpha,
        	    X->data, X->stride,
        	    Y->data, Y->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


/* SYR */

int  gsl_blas_ssyr (CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    gsl_matrix_float * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;
  if(X->size != row_dim) return GSL_EINVAL;

  gsl_blas_raw_ssyr(Uplo,
                    row_dim,
		    alpha,
        	    X->data, X->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


int  gsl_blas_dsyr (CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector * X,
                    gsl_matrix * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;
  if(X->size != row_dim) return GSL_EINVAL;

  gsl_blas_raw_dsyr(Uplo,
                    row_dim,
		    alpha,
        	    X->data, X->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


/* SPR */

int  gsl_blas_sspr (CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    float Ap[])
{
}


int  gsl_blas_dspr (CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector * X,
                    double Ap[])
{
}


/* SYR2 */

int  gsl_blas_ssyr2 (CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     gsl_matrix_float * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;
  if(X->size != col_dim) return GSL_EINVAL;
  if(Y->size != col_dim) return GSL_EINVAL;
  
  gsl_blas_raw_ssyr2(Uplo,
                     col_dim,
		     alpha,
		     X->data, X->stride,
		     Y->data, Y->stride,
		     A->data, tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_dsyr2 (CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     gsl_matrix * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->dim2 ;

  if(row_dim != col_dim) return GSL_EINVAL;
  if(X->size != col_dim) return GSL_EINVAL;
  if(Y->size != col_dim) return GSL_EINVAL;
  
  gsl_blas_raw_dsyr2(Uplo,
                     col_dim,
		     alpha,
		     X->data, X->stride,
		     Y->data, Y->stride,
		     A->data, tda
		     );
  return GSL_SUCCESS;
}


/* SPR2 */

int  gsl_blas_sspr2 (CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     gsl_vector_float * A)
{
}


int  gsl_blas_dspr2 (CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     double Ap[])
{
}


/* HEMV */

int  gsl_blas_chemv (CBLAS_UPLO Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
}


int  gsl_blas_zhemv (CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


/* HBMV */

int  gsl_blas_chbmv (CBLAS_UPLO Uplo,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
}


int  gsl_blas_zhbmv (CBLAS_UPLO Uplo,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


/* HPMV */

int  gsl_blas_chpmv (CBLAS_UPLO Uplo,
                     const gsl_complex_float * alpha,
		     const void * Ap,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
}


int  gsl_blas_zhpmv (CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const void * Ap,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
}


/* GERU */

int  gsl_blas_cgeru (const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A)
{
}


int  gsl_blas_zgeru (const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


/* GERC */

int  gsl_blas_cgerc (const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A)
{
}


int  gsl_blas_zgerc (const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


/* HER */

int  gsl_blas_cher (CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_complex_float * X,
                    gsl_matrix_complex_float * A)
{
}


int  gsl_blas_zher (CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    gsl_matrix_complex * A)
{
}


/* HPR */

int  gsl_blas_chpr (CBLAS_UPLO Uplo,
                    float alpha,
                    const gsl_vector_complex_float * X,
                    void * Ap)
{
}


int  gsl_blas_zhpr (CBLAS_UPLO Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    void * Ap)
{
}


/* HER2 */

int  gsl_blas_cher2 (CBLAS_UPLO Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A)
{
}


int  gsl_blas_zher2 (CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
}


/* HPR2 */

int  gsl_blas_chpr2 (CBLAS_UPLO Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     void * Ap)
{
}


int  gsl_blas_zhpr2 (CBLAS_UPLO Uplo,
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


/* GEMM */

int  gsl_blas_sgemm (CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C)
{
  if(A->size2 != B->size1) return GSL_EINVAL;
  if(A->size1 != C->size1) return GSL_EINVAL;
  if(B->size2 != C->size2) return GSL_EINVAL;
  
  gsl_blas_raw_sgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     alpha,
		     A->data, A->dim2,
		     B->data, B->dim2,
		     beta,
		     C->data, C->dim2
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_dgemm (CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C)
{
  if(A->size2 != B->size1) return GSL_EINVAL;
  if(A->size1 != C->size1) return GSL_EINVAL;
  if(B->size2 != C->size2) return GSL_EINVAL;
  
  gsl_blas_raw_dgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     alpha,
		     A->data, A->dim2,
		     B->data, B->dim2,
		     beta,
		     C->data, C->dim2
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_cgemm (CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
  if(A->size2 != B->size1) return GSL_EINVAL;
  if(A->size1 != C->size1) return GSL_EINVAL;
  if(B->size2 != C->size2) return GSL_EINVAL;
  
  gsl_blas_raw_cgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     (float *)alpha->dat,
		     A->data, A->dim2,
		     B->data, B->dim2,
		     (float *)beta->dat,
		     C->data, C->dim2
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_zgemm (CBLAS_TRANSPOSE TransA,
                     CBLAS_TRANSPOSE TransB,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
  if(A->size2 != B->size1) return GSL_EINVAL;
  if(A->size1 != C->size1) return GSL_EINVAL;
  if(B->size2 != C->size2) return GSL_EINVAL;
  
  gsl_blas_raw_zgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     (double *)alpha->dat,
		     A->data, A->dim2,
		     B->data, B->dim2,
		     (double *)beta->dat,
		     C->data, C->dim2
		     );
  return GSL_SUCCESS;
}


/* SYMM */

int  gsl_blas_ssymm (CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C)
{

}


int  gsl_blas_dsymm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C)
{
}


int  gsl_blas_csymm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
}

int  gsl_blas_zsymm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


/* SYRK */

int  gsl_blas_ssyrk (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     float beta,
                     gsl_matrix_float * C)
{
}


int  gsl_blas_dsyrk (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     double beta,
                     gsl_matrix * C)
{
}


int  gsl_blas_csyrk (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zsyrk (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


/* SYR2K */

int  gsl_blas_ssyr2k (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                      int K,
                      float alpha,
                      const gsl_matrix_float * A,
                      const gsl_matrix_float * B,
                      float beta,
                      gsl_matrix_float * C)
{
}


int  gsl_blas_dsyr2k (CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      double alpha,
                      const  gsl_matrix * A,
                      const  gsl_matrix * B,
                      double beta,
                      gsl_matrix * C)
{
}


int  gsl_blas_csyr2k (CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex_float * alpha,
                      const gsl_matrix_complex_float * A,
                      const gsl_matrix_complex_float * B,
                      const gsl_complex_float * beta,
                      gsl_matrix_complex_float * C)
{
}



int  gsl_blas_zsyr2k (CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      const gsl_complex * beta,
                      gsl_matrix_complex *C)
{
}


/* TRMM */

int  gsl_blas_strmm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B)
{
}


int  gsl_blas_dtrmm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B)
{
}


int  gsl_blas_ctrmm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     gsl_matrix_complex_float * B)
{
}


int  gsl_blas_ztrmm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


/* TRSM */

int  gsl_blas_strsm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B)
{
}


int  gsl_blas_dtrsm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B)
{
}


int  gsl_blas_ctrsm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     gsl_matrix_complex_float * B)
{
}


int  gsl_blas_ztrsm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                     CBLAS_DIAG Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


/* HEMM */

int  gsl_blas_chemm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zhemm (CBLAS_SIDE Side,
                     CBLAS_UPLO Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


/* HERK */

int  gsl_blas_cherk (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_complex_float * A,
                     float beta,
                     gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zherk (CBLAS_UPLO Uplo,
                     CBLAS_TRANSPOSE Trans,
                     int K,
                     double alpha,
                     const gsl_matrix_complex * A,
                     double beta,
                     gsl_matrix_complex * C)
{
}


/* HER2K */

int  gsl_blas_cher2k (CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex_float * alpha,
                      const gsl_matrix_complex_float * A,
                      const gsl_matrix_complex_float * B,
                      float beta,
                      gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zher2k (CBLAS_UPLO Uplo,
                      CBLAS_TRANSPOSE Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      double beta,
                      gsl_matrix_complex * C)
{
}
