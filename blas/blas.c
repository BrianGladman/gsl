/* blas/blas.c
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
/* GSL implementation of BLAS operations.
 * Conforms to gsl_blas interface.
 * Note that GSL native storage is row-major, so
 * we implement in terms of gsl_blas_raw.
 */
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
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
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
  *result = gsl_blas_raw_sdsdot(N, alpha, X->data, X->stride, Y->data, Y->stride);
  return status;
}

int gsl_blas_dsdot (const gsl_vector_float * X,
                    const gsl_vector_float * Y,
		    double * result
		    )
{
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
  *result = gsl_blas_raw_dsdot(N, X->data, X->stride, Y->data, Y->stride);
  return status;
}

int gsl_blas_sdot (const gsl_vector_float * X,
                   const gsl_vector_float * Y,
		   float * result
		   )
{
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
  *result = gsl_blas_raw_sdot(N, X->data, X->stride, Y->data, Y->stride);
  return status;
}

int gsl_blas_ddot (const gsl_vector * X,
                   const gsl_vector * Y,
		   double * result
		   )
{
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
  *result = gsl_blas_raw_ddot(N, X->data, X->stride, Y->data, Y->stride);
  return status;
}


int  gsl_blas_cdotu (const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_complex_float * dotu)
{
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
  gsl_blas_raw_cdotu(N, X->data, X->stride, Y->data, Y->stride, dotu->dat);
  return status;
}


int  gsl_blas_cdotc (const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_complex_float * dotc)
{
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
  gsl_blas_raw_cdotc(N, X->data, X->stride, Y->data, Y->stride, dotc->dat);
  return status;
}


int  gsl_blas_zdotu (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotu)
{
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
  gsl_blas_raw_zdotu(N, X->data, X->stride, Y->data, Y->stride, dotu->dat);
  return status;
}


int  gsl_blas_zdotc (const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_complex * dotc)
{
  CBLAS_INDEX_t N = GSL_MIN(X->size, Y->size);
  int status = GSL_SUCCESS;
  if(X->size != Y->size) status = GSL_EBADLEN;
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


CBLAS_INDEX_t gsl_blas_isamax (const gsl_vector_float * X)
{
  return gsl_blas_raw_isamax(X->size, X->data, X->stride);
}

CBLAS_INDEX_t gsl_blas_idamax (const gsl_vector * X)
{
  return gsl_blas_raw_idamax(X->size, X->data, X->stride);
}

CBLAS_INDEX_t gsl_blas_icamax (const gsl_vector_complex_float * X)
{
  return gsl_blas_raw_icamax(X->size, X->data, X->stride);
}

CBLAS_INDEX_t gsl_blas_izamax (const gsl_vector_complex * X)
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
    return GSL_EBADLEN;
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
    return GSL_EBADLEN;
}


int
gsl_blas_saxpy (float alpha, const gsl_vector_float * X, gsl_vector_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_saxpy(X->size, alpha, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_dswap (gsl_vector * X, gsl_vector * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_dswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_dcopy (const gsl_vector * X, gsl_vector * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_dcopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_daxpy (double alpha, const gsl_vector * X, gsl_vector * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_daxpy(X->size, alpha, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_cswap (gsl_vector_complex_float * X, gsl_vector_complex_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_cswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_ccopy (const gsl_vector_complex_float * X, gsl_vector_complex_float * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_ccopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
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
    return GSL_EBADLEN;
}


int
gsl_blas_zswap (gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_zswap(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_zcopy (const gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_zcopy(X->size, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_zaxpy (const gsl_complex * alpha, const gsl_vector_complex * X, gsl_vector_complex * Y)
{
  if(X->size == Y->size) {
    gsl_blas_raw_zaxpy(X->size, (double *)alpha->dat, X->data, X->stride, Y->data, Y->stride);
    return GSL_SUCCESS;
  }
  else
    return GSL_EBADLEN;
}


int
gsl_blas_srotg (float a[], float b[], float c[], float s[])
{
  gsl_blas_raw_srotg(a, b, c, s);
  return GSL_SUCCESS;
}


int
gsl_blas_srotmg (float d1[], float d2[], float b1[], float b2, float P[])
{
  gsl_blas_raw_srotmg(d1, d2, b1, b2, P);
  return GSL_SUCCESS;
}


int
gsl_blas_srot (gsl_vector_float * X, gsl_vector_float * Y, float c, float s)
{
  gsl_blas_raw_srot(GSL_MIN(X->size,Y->size), X->data, X->stride, Y->data, Y->stride, c, s);
  if(X->size != Y->size) {
    return GSL_EBADLEN;
  }
  else {
    return GSL_SUCCESS;
  }
}


int
gsl_blas_srotm (gsl_vector_float * X, gsl_vector_float * Y, const float P[])
{
  gsl_blas_raw_srotm(GSL_MIN(X->size,Y->size), X->data, X->stride, Y->data, Y->stride, P);
  if(X->size != Y->size) {
    return GSL_EBADLEN;
  }
  else {
    return GSL_SUCCESS;
  }
}


int
gsl_blas_drotg (double a[], double b[], double c[], double s[])
{
  gsl_blas_raw_drotg(a, b, c, s);
  return GSL_SUCCESS;
}


int
gsl_blas_drotmg (double d1[], double d2[], double b1[], double b2, double P[])
{
  gsl_blas_raw_drotmg(d1, d2, b1, b2, P);
  return GSL_SUCCESS;
}


int
gsl_blas_drot (gsl_vector * X, gsl_vector * Y, const double c, const double s)
{
  gsl_blas_raw_drot(GSL_MIN(X->size, Y->size), X->data, X->stride, Y->data, Y->stride, c, s);
  if(X->size != Y->size) {
    return GSL_EBADLEN;
  }
  else {
    return GSL_SUCCESS;
  }
}


int
gsl_blas_drotm (gsl_vector * X, gsl_vector * Y, const double P[])
{
  gsl_blas_raw_drotm(GSL_MIN(X->size, Y->size), X->data, X->stride, Y->data, Y->stride, P);
  if(X->size != Y->size) {
    return GSL_EBADLEN;
  }
  else {
    return GSL_SUCCESS;
  }
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

int  gsl_blas_sgemv (CBLAS_TRANSPOSE_t TransA,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}


int  gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}


int  gsl_blas_cgemv (CBLAS_TRANSPOSE_t TransA,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}


int  gsl_blas_zgemv (CBLAS_TRANSPOSE_t TransA,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}


/* GBMV */

int  gsl_blas_sgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
		     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}

int  gsl_blas_dgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}


int  gsl_blas_cgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}


int  gsl_blas_zgbmv (CBLAS_TRANSPOSE_t TransA,
                     int KL, int KU,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if((TransA == CblasNoTrans && col_dim == X->size && row_dim == Y->size) ||
     (TransA == CblasTrans   && row_dim == X->size && col_dim == Y->size)
     ) {
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
    return GSL_EBADLEN;
}


/* TRMV */

int  gsl_blas_strmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_strmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_dtrmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_dtrmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ctrmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ctrmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ztrmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ztrmv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


/* TBMV */

int  gsl_blas_stbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_stbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_dtbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_dtbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ctbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ctbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ztbmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ztbmv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


/* TPMV */

int  gsl_blas_stpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const float Ap[],
                     gsl_vector_float * X)
{
  gsl_blas_raw_stpmv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_dtpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const double * Ap,
                     gsl_vector * X)
{
  gsl_blas_raw_dtpmv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_ctpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void * Ap,
                     gsl_vector_complex_float * X)
{
  gsl_blas_raw_ctpmv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_ztpmv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void * Ap,
                     gsl_vector_complex * X)
{
  gsl_blas_raw_ztpmv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


/* TRSV */

int  gsl_blas_strsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_strsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_dtrsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_dtrsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ctrsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ctrsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ztrsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex *X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ztrsv(Uplo, TransA, Diag,
                       row_dim,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


/* TBSV */

int  gsl_blas_stbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_float * A,
                     gsl_vector_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_stbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_dtbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix * A,
                     gsl_vector * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_dtbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ctbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex_float * A,
                     gsl_vector_complex_float * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ctbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_ztbsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     int K,
                     const gsl_matrix_complex * A,
                     gsl_vector_complex * X)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ztbsv(Uplo, TransA, Diag,
                       row_dim,
		       K,
		       A->data, tda,
		       X->data, X->stride
		       );
    return GSL_SUCCESS;
  }
}


/* TPSV */

int  gsl_blas_stpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const float Ap[],
                     gsl_vector_float * X)
{
  gsl_blas_raw_stpsv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_dtpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const double Ap[],
                     gsl_vector * X)
{
  gsl_blas_raw_dtpsv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_ctpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void * Ap,
                     gsl_vector_complex_float * X)
{
  gsl_blas_raw_ctpsv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_ztpsv (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                     const void *Ap,
                     gsl_vector_complex * X)
{
  gsl_blas_raw_ztpsv(Uplo, TransA, Diag, X->size, Ap, X->data, X->stride);
  return GSL_SUCCESS;
}


/* SYMV */

int  gsl_blas_ssymv (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
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
}


int  gsl_blas_dsymv (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
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
}


/* SBMV */

int  gsl_blas_ssbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
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
}


int  gsl_blas_dsbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(col_dim != X->size) {
    return GSL_EBADLEN;
  }
  else {
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
}


/* SPMV */

int  gsl_blas_sspmv (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const float * Ap,
                     const gsl_vector_float * X,
                     float beta,
                     gsl_vector_float * Y)
{
  gsl_blas_raw_sspmv(Uplo, X->size, alpha, Ap, X->data, X->stride, beta, Y->data, Y->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_dspmv (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const double * Ap,
                     const gsl_vector * X,
                     double beta,
                     gsl_vector * Y)
{
  gsl_blas_raw_dspmv(Uplo, X->size, alpha, Ap, X->data, X->stride, beta, Y->data, Y->stride);
  return GSL_SUCCESS;
}


/* GER */

int  gsl_blas_sger (float alpha,
                    const gsl_vector_float * X,
                    const gsl_vector_float * Y,
                    gsl_matrix_float * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->tda ;

  /* FIXME: check this */
  if(X->size != row_dim || Y->size != col_dim) return GSL_EBADLEN;

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
  size_t tda = A->tda ;

  /* FIXME: check this */
  if(X->size != row_dim || Y->size != col_dim) return GSL_EBADLEN;

  gsl_blas_raw_dger(row_dim, col_dim,
		    alpha,
        	    X->data, X->stride,
        	    Y->data, Y->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


/* SYR */

int  gsl_blas_ssyr (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    gsl_matrix_float * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) {
    return GSL_ENOTSQR;
  }
  else if(X->size != row_dim) {
    return GSL_EBADLEN;
  }
  else {
    gsl_blas_raw_ssyr(Uplo,
                      row_dim,
	  	      alpha,
        	      X->data, X->stride,
        	      A->data, tda
        	      );
    return GSL_SUCCESS;
  }
}


int  gsl_blas_dsyr (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector * X,
                    gsl_matrix * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(X->size != row_dim) return GSL_EBADLEN;

  gsl_blas_raw_dsyr(Uplo,
                    row_dim,
		    alpha,
        	    X->data, X->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


/* SPR */

int  gsl_blas_sspr (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_float * X,
                    float * Ap)
{
  gsl_blas_raw_sspr(Uplo, X->size, alpha, X->data, X->stride, Ap);
  return GSL_SUCCESS;
}


int  gsl_blas_dspr (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector * X,
                    double * Ap)
{
  gsl_blas_raw_dspr(Uplo, X->size, alpha, X->data, X->stride, Ap);
  return GSL_SUCCESS;
}


/* SYR2 */

int  gsl_blas_ssyr2 (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
                     gsl_matrix_float * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(X->size != col_dim) return GSL_EBADLEN;
  if(Y->size != col_dim) return GSL_EBADLEN;
  
  gsl_blas_raw_ssyr2(Uplo,
                     col_dim,
		     alpha,
		     X->data, X->stride,
		     Y->data, Y->stride,
		     A->data, tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_dsyr2 (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     gsl_matrix * A)
{
  size_t  row_dim = A->size1;
  size_t  col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(X->size != col_dim) return GSL_EBADLEN;
  if(Y->size != col_dim) return GSL_EBADLEN;
  
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

int  gsl_blas_sspr2 (CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_vector_float * X,
                     const gsl_vector_float * Y,
		     float * Ap)
{
  gsl_blas_raw_sspr2(Uplo, X->size, alpha, X->data, X->stride, Y->data, Y->stride, Ap);
  return GSL_SUCCESS;
}


int  gsl_blas_dspr2 (CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_vector * X,
                     const gsl_vector * Y,
                     double * Ap)
{
  gsl_blas_raw_dspr2(Uplo, X->size, alpha, X->data, X->stride, Y->data, Y->stride, Ap);
  return GSL_SUCCESS;
}


/* HEMV */

int  gsl_blas_chemv (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(col_dim != X->size) return GSL_EBADLEN;

  gsl_blas_raw_chemv(Uplo,
 		     col_dim,
        	     alpha->dat,
        	     A->data, tda,
        	     X->data, X->stride,
        	     beta->dat,
        	     Y->data, Y->stride
        	     );
  return GSL_SUCCESS;
}


int  gsl_blas_zhemv (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(col_dim != X->size) return GSL_EBADLEN;

  gsl_blas_raw_zhemv(Uplo,
  		     col_dim,
        	     alpha->dat,
        	     A->data, tda,
        	     X->data, X->stride,
        	     beta->dat,
        	     Y->data, Y->stride
        	     );
  return GSL_SUCCESS;
}


/* HBMV */

int  gsl_blas_chbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(col_dim != X->size) return GSL_EBADLEN;

  gsl_blas_raw_chbmv(Uplo,
  		     K,
  		     col_dim,
        	     alpha->dat,
        	     A->data, tda,
        	     X->data, X->stride,
        	     beta->dat,
        	     Y->data, Y->stride
        	     );
  return GSL_SUCCESS;
}


int  gsl_blas_zhbmv (CBLAS_UPLO_t Uplo,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
  size_t row_dim  = A->size1;
  size_t col_dim  = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(col_dim != X->size) return GSL_EBADLEN;

  gsl_blas_raw_zhbmv(Uplo,
  		     K,
  		     col_dim,
        	     alpha->dat,
        	     A->data, tda,
        	     X->data, X->stride,
        	     beta->dat,
        	     Y->data, Y->stride
        	     );
  return GSL_SUCCESS;
}


/* HPMV */

int  gsl_blas_chpmv (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
		     const void * Ap,
                     const gsl_vector_complex_float * X,
                     const gsl_complex_float * beta,
                     gsl_vector_complex_float * Y)
{
  gsl_blas_raw_chpmv(Uplo, X->size, alpha->dat, Ap, X->data, X->stride, beta->dat, Y->data, Y->stride);
  return GSL_SUCCESS;
}


int  gsl_blas_zhpmv (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const void * Ap,
                     const gsl_vector_complex * X,
                     const gsl_complex * beta,
                     gsl_vector_complex * Y)
{
  gsl_blas_raw_zhpmv(Uplo, X->size, alpha->dat, Ap, X->data, X->stride, beta->dat, Y->data, Y->stride);
  return GSL_SUCCESS;
}


/* GERU */

int  gsl_blas_cgeru (const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A)
{
  gsl_blas_raw_cgeru (A->size1, A->size2, alpha->dat,
                      X->data, X->stride,
                      Y->data, Y->stride,
                      A->data, A->tda);
  return GSL_SUCCESS;
}


int  gsl_blas_zgeru (const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
  gsl_blas_raw_zgeru (A->size1, A->size2, alpha->dat,
                      X->data, X->stride,
                      Y->data, Y->stride,
                      A->data, A->tda);
  return GSL_SUCCESS;
}


/* GERC */

int  gsl_blas_cgerc (const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A)
{
  gsl_blas_raw_cgerc (A->size1, A->size2, alpha->dat,
                      X->data, X->stride,
                      Y->data, Y->stride,
                      A->data, A->tda);
  return GSL_SUCCESS;
}


int  gsl_blas_zgerc (const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
  gsl_blas_raw_zgerc (A->size1, A->size2, alpha->dat,
                      X->data, X->stride,
                      Y->data, Y->stride,
                      A->data, A->tda);
  return GSL_SUCCESS;
}


/* HER */

int  gsl_blas_cher (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_complex_float * X,
                    gsl_matrix_complex_float * A)
{
  size_t row_dim = A->size1;
  size_t col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(X->size != row_dim) return GSL_EBADLEN;

  gsl_blas_raw_cher(Uplo,
                    row_dim,
		    alpha,
        	    X->data, X->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


int  gsl_blas_zher (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    gsl_matrix_complex * A)
{
  size_t row_dim = A->size1;
  size_t col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(X->size != row_dim) return GSL_EBADLEN;

  gsl_blas_raw_zher(Uplo,
                    row_dim,
		    alpha,
        	    X->data, X->stride,
        	    A->data, tda
        	    );
  return GSL_SUCCESS;
}


/* HPR */

int  gsl_blas_chpr (CBLAS_UPLO_t Uplo,
                    float alpha,
                    const gsl_vector_complex_float * X,
                    void * Ap)
{
  gsl_blas_raw_chpr(Uplo, X->size, alpha, X->data, X->stride, Ap);
  return GSL_SUCCESS;
}


int  gsl_blas_zhpr (CBLAS_UPLO_t Uplo,
                    double alpha,
                    const gsl_vector_complex * X,
                    void * Ap)
{
  gsl_blas_raw_zhpr(Uplo, X->size, alpha, X->data, X->stride, Ap);
  return GSL_SUCCESS;
}


/* HER2 */

int  gsl_blas_cher2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     gsl_matrix_complex_float * A)
{
  size_t row_dim = A->size1;
  size_t col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(X->size != col_dim) return GSL_EBADLEN;
  if(Y->size != col_dim) return GSL_EBADLEN;
  
  gsl_blas_raw_cher2(Uplo,
                     col_dim,
		     alpha->dat,
		     X->data, X->stride,
		     Y->data, Y->stride,
		     A->data, tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_zher2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     gsl_matrix_complex * A)
{
  size_t row_dim = A->size1;
  size_t col_dim = A->size2;
  size_t tda = A->tda ;

  if(row_dim != col_dim) return GSL_ENOTSQR;
  if(X->size != col_dim) return GSL_EBADLEN;
  if(Y->size != col_dim) return GSL_EBADLEN;
  
  gsl_blas_raw_zher2(Uplo,
                     col_dim,
		     alpha->dat,
		     X->data, X->stride,
		     Y->data, Y->stride,
		     A->data, tda
		     );
  return GSL_SUCCESS;
}


/* HPR2 */

int  gsl_blas_chpr2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_vector_complex_float * X,
                     const gsl_vector_complex_float * Y,
                     void * Ap)
{
  gsl_blas_raw_chpr2(Uplo, X->size, alpha->dat, X->data, X->stride, Y->data, Y->stride, Ap);
  return GSL_SUCCESS;
}


int  gsl_blas_zhpr2 (CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_vector_complex * X,
                     const gsl_vector_complex * Y,
                     void * Ap)
{
  gsl_blas_raw_zhpr2(Uplo, X->size, alpha->dat, X->data, X->stride, Y->data, Y->stride, Ap);
  return GSL_SUCCESS;
}


/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */


/* GEMM */

int  gsl_blas_sgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;
  
  gsl_blas_raw_sgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     alpha,
		     A->data, A->tda,
		     B->data, B->tda,
		     beta,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_dgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;
  
  gsl_blas_raw_dgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     alpha,
		     A->data, A->tda,
		     B->data, B->tda,
		     beta,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_cgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;
  
  gsl_blas_raw_cgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     (float *)alpha->dat,
		     A->data, A->tda,
		     B->data, B->tda,
		     (float *)beta->dat,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_zgemm (CBLAS_TRANSPOSE_t TransA,
                     CBLAS_TRANSPOSE_t TransB,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;
  
  gsl_blas_raw_zgemm(TransA, TransB,
                     C->size1, C->size2,
		     K,
		     (double *)alpha->dat,
		     A->data, A->tda,
		     B->data, B->tda,
		     (double *)beta->dat,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}


/* SYMM */

int  gsl_blas_ssymm (CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo,
                     float alpha,
                     const gsl_matrix_float * A,
                     const gsl_matrix_float * B,
                     float beta,
                     gsl_matrix_float * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;

  gsl_blas_raw_ssymm(Side, Uplo,
                     C->size1, C->size2,
		     alpha,
		     A->data, A->tda,
		     B->data, B->tda,
		     beta,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_dsymm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     double alpha,
                     const gsl_matrix * A,
                     const gsl_matrix * B,
                     double beta,
                     gsl_matrix * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;

  gsl_blas_raw_dsymm(Side, Uplo,
                     C->size1, C->size2,
		     alpha,
		     A->data, A->tda,
		     B->data, B->tda,
		     beta,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}


int  gsl_blas_csymm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;

  gsl_blas_raw_csymm(Side, Uplo,
                     C->size1, C->size2,
		     alpha->dat,
		     A->data, A->tda,
		     B->data, B->tda,
		     beta->dat,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}

int  gsl_blas_zsymm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
  if(A->size2 != B->size1) return GSL_EBADLEN;
  if(A->size1 != C->size1) return GSL_EBADLEN;
  if(B->size2 != C->size2) return GSL_EBADLEN;

  gsl_blas_raw_zsymm(Side, Uplo,
                     C->size1, C->size2,
		     alpha->dat,
		     A->data, A->tda,
		     B->data, B->tda,
		     beta->dat,
		     C->data, C->tda
		     );
  return GSL_SUCCESS;
}


#if 0

/* SYRK */

int  gsl_blas_ssyrk (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_float * A,
                     float beta,
                     gsl_matrix_float * C)
{
}


int  gsl_blas_dsyrk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     double alpha,
                     const gsl_matrix * A,
                     double beta,
                     gsl_matrix * C)
{
}


int  gsl_blas_csyrk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zsyrk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


/* SYR2K */

int  gsl_blas_ssyr2k (CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t Trans,
                      int K,
                      float alpha,
                      const gsl_matrix_float * A,
                      const gsl_matrix_float * B,
                      float beta,
                      gsl_matrix_float * C)
{
}


int  gsl_blas_dsyr2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      double alpha,
                      const  gsl_matrix * A,
                      const  gsl_matrix * B,
                      double beta,
                      gsl_matrix * C)
{
}


int  gsl_blas_csyr2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex_float * alpha,
                      const gsl_matrix_complex_float * A,
                      const gsl_matrix_complex_float * B,
                      const gsl_complex_float * beta,
                      gsl_matrix_complex_float * C)
{
}



int  gsl_blas_zsyr2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      const gsl_complex * beta,
                      gsl_matrix_complex *C)
{
}


/* TRMM */

int  gsl_blas_strmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B)
{
}


int  gsl_blas_dtrmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B)
{
}


int  gsl_blas_ctrmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     gsl_matrix_complex_float * B)
{
}


int  gsl_blas_ztrmm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


/* TRSM */

int  gsl_blas_strsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     float alpha,
                     const gsl_matrix_float * A,
                     gsl_matrix_float * B)
{
}


int  gsl_blas_dtrsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     double alpha,
                     const gsl_matrix * A,
                     gsl_matrix * B)
{
}


int  gsl_blas_ctrsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     gsl_matrix_complex_float * B)
{
}


int  gsl_blas_ztrsm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA,
                     CBLAS_DIAG_t Diag,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     gsl_matrix_complex * B)
{
}


/* HEMM */

int  gsl_blas_chemm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex_float * alpha,
                     const gsl_matrix_complex_float * A,
                     const gsl_matrix_complex_float * B,
                     const gsl_complex_float * beta,
                     gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zhemm (CBLAS_SIDE_t Side,
                     CBLAS_UPLO_t Uplo,
                     const gsl_complex * alpha,
                     const gsl_matrix_complex * A,
                     const gsl_matrix_complex * B,
                     const gsl_complex * beta,
                     gsl_matrix_complex * C)
{
}


/* HERK */

int  gsl_blas_cherk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     float alpha,
                     const gsl_matrix_complex_float * A,
                     float beta,
                     gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zherk (CBLAS_UPLO_t Uplo,
                     CBLAS_TRANSPOSE_t Trans,
                     int K,
                     double alpha,
                     const gsl_matrix_complex * A,
                     double beta,
                     gsl_matrix_complex * C)
{
}


/* HER2K */

int  gsl_blas_cher2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex_float * alpha,
                      const gsl_matrix_complex_float * A,
                      const gsl_matrix_complex_float * B,
                      float beta,
                      gsl_matrix_complex_float * C)
{
}


int  gsl_blas_zher2k (CBLAS_UPLO_t Uplo,
                      CBLAS_TRANSPOSE_t Trans,
                      int K,
                      const gsl_complex * alpha,
                      const gsl_matrix_complex * A,
                      const gsl_matrix_complex * B,
                      double beta,
                      gsl_matrix_complex * C)
{
}

#endif /* 0 */
