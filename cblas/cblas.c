/* blas/blas_raw_native.c
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
 * Author:  G. Jungman  and  B. Gough
 * RCS:     $Id$
 */

/* 
 * Conforms to cblas_ interface.
 */

#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#define OFFSET(N, incX) ((incX) > 0 ?  0 : ((N) - 1) * (-(incX)))
#define BLAS_ERROR(x)  abort(); 

#define CONJUGATE(x) ((x) == CblasConjTrans)
#define TRANSPOSE(x) ((x) == CblasTrans || (x) == CblasConjTrans)
#define UPPER(x) ((x) == CblasUpper)
#define LOWER(x) ((x) == CblasLower)

/* Handling of packed complex types... */

#define REAL(a,i) (((BASE *) a)[2*(i)])
#define IMAG(a,i) (((BASE *) a)[2*(i)+1])

#define REAL0(a) (((BASE *)a)[0])
#define IMAG0(a) (((BASE *)a)[1])

#define GB(KU,KL,lda,i,j) ((KU+1+(i-j))*lda + j)

#define TRCOUNT(N,i) ((((i)+1)*(2*(N)-(i)))/2)

#define TBUP(N,i,j) 
#define TBLO(N,i,j) 

#define TPUP(N,i,j) (TRCOUNT(N,(i)-1)+(j)-(i))
#define TPLO(N,i,j) (((i)*((i)+1))/2 + (j))

/* ===========================================================================
 * Level 1 BLAS functions
 * ===========================================================================
 */

float
cblas_sdsdot (const int N, const float alpha, const float *X, const int incX,
	      const float *Y, const int incY)
{
#define INIT_VAL  alpha
#define ACC_TYPE  double
#define BASE float
#include "source_dot_r.h"
#undef ACC_TYPE
#undef BASE
#undef INIT_VAL
}

double
cblas_dsdot (const int N, const float *X, const int incX, const float *Y,
	     const int incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  double
#define BASE float
#include "source_dot_r.h"
#undef ACC_TYPE
#undef BASE
#undef INIT_VAL
}

float
cblas_sdot (const int N, const float *X, const int incX, const float *Y,
	    const int incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  float
#define BASE float
#include "source_dot_r.h"
#undef ACC_TYPE
#undef BASE
#undef INIT_VAL
}

double
cblas_ddot (const int N, const double *X, const int incX, const double *Y,
	    const int incY)
{
#define INIT_VAL  0.0
#define ACC_TYPE  double
#define BASE double
#include "source_dot_r.h"
#undef ACC_TYPE
#undef BASE
#undef INIT_VAL
}


void
cblas_cdotu_sub (const int N, const void *X, const int incX, const void *Y,
	     const int incY, void *result)
{
#define BASE float
#define CONJ_SIGN 1.0
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE
}

void
cblas_cdotc_sub (const int N, const void *X, const int incX, const void *Y,
	     const int incY, void *result)
{
#define BASE float
#define CONJ_SIGN (-1.0)
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE
}

void
cblas_zdotu_sub (const int N, const void *X, const int incX, const void *Y,
	     const int incY, void *result)
{
#define BASE double
#define CONJ_SIGN 1.0
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE
}

void
cblas_zdotc_sub (const int N, const void *X, const int incX, const void *Y,
	     const int incY, void *result)
{
#define BASE double
#define CONJ_SIGN (-1.0)
#include "source_dot_c.h"
#undef CONJ_SIGN
#undef BASE
}


float
cblas_snrm2 (const int N, const float *X, const int incX)
{
#define BASE float
#include "source_nrm2_r.h"
#undef BASE
}

double
cblas_dnrm2 (const int N, const double *X, const int incX)
{
#define BASE double
#include "source_nrm2_r.h"
#undef BASE
}

float
cblas_scnrm2 (const int N, const void *X, const int incX)
{
#define BASE float
#include "source_nrm2_c.h"
#undef BASE
}

double
cblas_dznrm2 (const int N, const void *X, const int incX)
{
#define BASE double
#include "source_nrm2_c.h"
#undef BASE
}


float
cblas_sasum (const int N, const float *X, const int incX)
{
#define BASE float
#include "source_asum_r.h"
#undef BASE
}

double
cblas_dasum (const int N, const double *X, const int incX)
{
#define BASE double
#include "source_asum_r.h"
#undef BASE
}

float
cblas_scasum (const int N, const void *X, const int incX)
{
#define BASE float
#include "source_asum_c.h"
#undef BASE
}

double
cblas_dzasum (const int N, const void *X, const int incX)
{
#define BASE double
#include "source_asum_c.h"
#undef BASE
}


CBLAS_INDEX
cblas_isamax (const int N, const float *X, const int incX)
{
#define BASE float
#include "source_iamax_r.h"
#undef BASE
}

CBLAS_INDEX
cblas_idamax (const int N, const double *X, const int incX)
{
#define BASE double
#include "source_iamax_r.h"
#undef BASE
}

CBLAS_INDEX
cblas_icamax (const int N, const void *X, const int incX)
{
#define BASE float
#include "source_iamax_c.h"
#undef BASE
}

CBLAS_INDEX
cblas_izamax (const int N, const void *X, const int incX)
{
#define BASE double
#include "source_iamax_c.h"
#undef BASE
}


void
cblas_sswap (const int N, float *X, const int incX, float *Y, const int incY)
{
#define BASE float
#include "source_swap_r.h"
#undef BASE
}

void
cblas_dswap (const int N, double *X, const int incX, double *Y,
	     const int incY)
{
#define BASE double
#include "source_swap_r.h"
#undef BASE
}

void
cblas_cswap (const int N, void *X, const int incX, void *Y, const int incY)
{
#define BASE float
#include "source_swap_c.h"
#undef BASE
}

void
cblas_zswap (const int N, void *X, const int incX, void *Y, const int incY)
{
#define BASE double
#include "source_swap_c.h"
#undef BASE
}


void
cblas_scopy (const int N, const float *X, const int incX, float *Y,
	     const int incY)
{
#define BASE float
#include "source_copy_r.h"
#undef BASE
}

void
cblas_dcopy (const int N, const double *X, const int incX, double *Y,
	     const int incY)
{
#define BASE double
#include "source_copy_r.h"
#undef BASE
}

void
cblas_ccopy (const int N, const void *X, const int incX, void *Y,
	     const int incY)
{
#define BASE float
#include "source_copy_c.h"
#undef BASE
}

void
cblas_zcopy (const int N, const void *X, const int incX, void *Y,
	     const int incY)
{
#define BASE double
#include "source_copy_c.h"
#undef BASE
}


void
cblas_saxpy (const int N, const float alpha, const float *X, const int incX,
	     float *Y, const int incY)
{
#define BASE float
#include "source_axpy_r.h"
#undef BASE
}

void
cblas_daxpy (const int N, const double alpha, const double *X, const int incX,
	     double *Y, const int incY)
{
#define BASE double
#include "source_axpy_r.h"
#undef BASE
}

void
cblas_caxpy (const int N, const void *alpha, const void *X, const int incX,
	     void *Y, const int incY)
{
#define BASE float
#include "source_axpy_c.h"
#undef BASE
}

void
cblas_zaxpy (const int N, const void *alpha, const void *X, const int incX,
	     void *Y, const int incY)
{
#define BASE double
#include "source_axpy_c.h"
#undef BASE
}


void
cblas_srotg (float *a, float *b, float *c, float *s)
{
#define BASE float
#include "source_rotg.h"
#undef BASE
}

void
cblas_drotg (double *a, double *b, double *c, double *s)
{
#define BASE double
#include "source_rotg.h"
#undef BASE
}

void
cblas_srotmg (float *d1, float *d2, float *b1, const float b2, float *P)
{
#define BASE float
#include "source_rotmg.h"
#undef BASE
}

void
cblas_drotmg (double *d1, double *d2, double *b1, double b2, double *P)
{
#define BASE double
#include "source_rotmg.h"
#undef BASE
}


void
cblas_srot (const int N, float *X, const int incX, float *Y, const int incY,
	    const float c, const float s)
{
#define BASE float
#include "source_rot.h"
#undef BASE
}

void
cblas_drot (const int N, double *X, const int incX, double *Y, const int incY,
	    const double c, const double s)
{
#define BASE double
#include "source_rot.h"
#undef BASE
}


void
cblas_srotm (const int N, float *X, const int incX, float *Y, const int incY,
	     const float *P)
{
#define BASE float
#include "source_rotm.h"
#undef BASE
}

void
cblas_drotm (const int N, double *X, const int incX, double *Y,
	     const int incY, const double *P)
{
#define BASE double
#include "source_rotm.h"
#undef BASE
}


void
cblas_sscal (const int N, const float alpha, float *X, const int incX)
{
#define BASE float
#include "source_scal_r.h"
#undef BASE
}

void
cblas_dscal (const int N, const double alpha, double *X, const int incX)
{
#define BASE double
#include "source_scal_r.h"
#undef BASE
}

void
cblas_cscal (const int N, const void *alpha, void *X, const int incX)
{
#define BASE float
#include "source_scal_c.h"
#undef BASE
}

void
cblas_zscal (const int N, const void *alpha, void *X, const int incX)
{
#define BASE double
#include "source_scal_c.h"
#undef BASE
}

void
cblas_csscal (const int N, const float alpha, void *X, const int incX)
{
#define BASE float
#include "source_scal_c_s.h"
#undef BASE
}

void
cblas_zdscal (const int N, const double alpha, void *X, const int incX)
{
#define BASE double
#include "source_scal_c_s.h"
#undef BASE
}


/*
 * ===========================================================================
 * Level 2 BLAS
 * ===========================================================================
 */


/* GEMV */

void
cblas_sgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const float alpha, const float *A,
	     const int lda, const float *X, const int incX, const float beta,
	     float *Y, const int incY)
{
#define BASE float
#include "source_gemv_r.h"
#undef BASE
}

void
cblas_dgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const double alpha, const double *A,
	     const int lda, const double *X, const int incX,
	     const double beta, double *Y, const int incY)
{
#define BASE double
#include "source_gemv_r.h"
#undef BASE
}

void
cblas_cgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const void *alpha, const void *A,
	     const int lda, const void *X, const int incX, const void *beta,
	     void *Y, const int incY)
{
#define BASE float
#include "source_gemv_c.h"
#undef BASE
}

void
cblas_zgemv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const void *alpha, const void *A,
	     const int lda, const void *X, const int incX, const void *beta,
	     void *Y, const int incY)
{
#define BASE double
#include "source_gemv_c.h"
#undef BASE
}

/* GBMV */

void
cblas_sgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const int KL, const int KU,
	     const float alpha, const float *A, const int lda, const float *X,
	     const int incX, const float beta, float *Y, const int incY)
{
#define BASE float
#include "source_gbmv_r.h"
#undef BASE
}

void
cblas_dgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const int KL, const int KU,
	     const double alpha, const double *A, const int lda,
	     const double *X, const int incX, const double beta, double *Y,
	     const int incY)
{
#define BASE double
#include "source_gbmv_r.h"
#undef BASE
}

void
cblas_cgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const int KL, const int KU,
	     const void *alpha, const void *A, const int lda, const void *X,
	     const int incX, const void *beta, void *Y, const int incY)
{
#define BASE float
#include "source_gbmv_c.h"
#undef BASE
}

void
cblas_zgbmv (const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA,
	     const int M, const int N, const int KL, const int KU,
	     const void *alpha, const void *A, const int lda, const void *X,
	     const int incX, const void *beta, void *Y, const int incY)
{
#define BASE double
#include "source_gbmv_c.h"
#undef BASE
}

/* HEMV */

void
cblas_chemv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *A, const int lda,
	     const void *X, const int incX, const void *beta, void *Y,
	     const int incY)
{
#define BASE float
#include "source_hemv.h"
#undef BASE
}

void
cblas_zhemv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *A, const int lda,
	     const void *X, const int incX, const void *beta, void *Y,
	     const int incY)
{
#define BASE double
#include "source_hemv.h"
#undef BASE
}


/* HBMV */

void
cblas_chbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const int K, const void *alpha, const void *A,
	     const int lda, const void *X, const int incX, const void *beta,
	     void *Y, const int incY)
{
#define BASE float
#include "source_hbmv.h"
#undef BASE
}

void
cblas_zhbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const int K, const void *alpha, const void *A,
	     const int lda, const void *X, const int incX, const void *beta,
	     void *Y, const int incY)
{
#define BASE double
#include "source_hbmv.h"
#undef BASE
}


/* HPMV */

void
cblas_chpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *Ap, const void *X,
	     const int incX, const void *beta, void *Y, const int incY)
{
#define BASE float
#include "source_hpmv.h"
#undef BASE
}

void
cblas_zhpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *Ap, const void *X,
	     const int incX, const void *beta, void *Y, const int incY)
{
#define BASE double
#include "source_hpmv.h"
#undef BASE
}

#ifdef 0
/* SYMV */

void
cblas_ssymv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const float alpha, const float *A, const int lda,
	     const float *X, const int incX, const float beta, float *Y,
	     const int incY)
{
#define BASE float
#include "source_symv.h"
#undef BASE
}

void
cblas_dsymv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const double alpha, const double *A, const int lda,
	     const double *X, const int incX, const double beta, double *Y,
	     const int incY)
{
#define BASE double
#include "source_symv.h"
#undef BASE
}


/* SBMV */

void
cblas_ssbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const int K, const float alpha, const float *A,
	     const int lda, const float *X, const int incX, const float beta,
	     float *Y, const int incY)
{
#define BASE float
#include "source_sbmv.h"
#undef BASE
}

void
cblas_dsbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const int K, const double alpha, const double *A,
	     const int lda, const double *X, const int incX,
	     const double beta, double *Y, const int incY)
{
#define BASE double
#include "source_sbmv.h"
#undef BASE
}

/* SPMV */

void
cblas_sspmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const float alpha, const float *Ap, const float *X,
	     const int incX, const float beta, float *Y, const int incY)
{
#define BASE float
#include "source_spmv.h"
#undef BASE
}

void
cblas_dspmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const double alpha, const double *Ap,
	     const double *X, const int incX, const double beta, double *Y,
	     const int incY)
{
#define BASE double
#include "source_spmv.h"
#undef BASE
}


/* TRMV */

void
cblas_strmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const float *A, const int lda, float *X,
	     const int incX)
{
#define BASE float
#include "source_trmv_r.h"
#undef BASE
}

void
cblas_dtrmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const double *A, const int lda, double *X,
	     const int incX)
{
#define BASE double
#include "source_trmv_r.h"
#undef BASE
}

void
cblas_ctrmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE float
#include "source_trmv_c.h"
#undef BASE
}

void
cblas_ztrmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE double
#include "source_trmv_c.h"
#undef BASE
}


/* TBMV */

void
cblas_stbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const float *A, const int lda,
	     float *X, const int incX)
{
#define BASE float
#include "source_tbmv_r.h"
#undef BASE
}

void
cblas_dtbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const double *A, const int lda,
	     double *X, const int incX)
{
#define BASE double
#include "source_tbmv_r.h"
#undef BASE
}

void
cblas_ctbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE float
#include "source_tbmv_c.h"
#undef BASE
}

void
cblas_ztbmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE double
#include "source_tbmv_c.h"
#undef BASE
}


/* TPMV */

void
cblas_stpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const float *Ap, float *X, const int incX)
{
#define BASE float
#include "source_tpmv_r.h"
#undef BASE
}

void
cblas_dtpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const double *Ap, double *X, const int incX)
{
#define BASE double
#include "source_tpmv_r.h"
#undef BASE
}

void
cblas_ctpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *Ap, void *X, const int incX)
{
#define BASE float
#include "source_tpmv_c.h"
#undef BASE
}

void
cblas_ztpmv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *Ap, void *X, const int incX)
{
#define BASE double
#include "source_tpmv_c.h"
#undef BASE
}


/* TRSV */

void
cblas_strsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const float *A, const int lda, float *X,
	     const int incX)
{
#define BASE float
#include "source_trsv_r.h"
#undef BASE
}

void
cblas_dtrsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const double *A, const int lda, double *X,
	     const int incX)
{
#define BASE double
#include "source_trsv_r.h"
#undef BASE
}

void
cblas_ctrsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE float
#include "source_trsv_c.h"
#undef BASE
}

void
cblas_ztrsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE double
#include "source_trsv_c.h"
#undef BASE
}


/* TBSV */

void
cblas_stbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const float *A, const int lda,
	     float *X, const int incX)
{
#define BASE float
#include "source_tbsv_r.h"
#undef BASE
}

void
cblas_dtbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const double *A, const int lda,
	     double *X, const int incX)
{
#define BASE double
#include "source_tbsv_r.h"
#undef BASE
}

void
cblas_ctbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE float
#include "source_tbsv_c.h"
#undef BASE
}

void
cblas_ztbsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const int K, const void *A, const int lda, void *X,
	     const int incX)
{
#define BASE double
#include "source_tbsv_c.h"
#undef BASE
}


/* TPSV */

void
cblas_stpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const float *Ap, float *X, const int incX)
{
#define BASE float
#include "source_tpsv_r.h"
#undef BASE
}

void
cblas_dtpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const double *Ap, double *X, const int incX)
{
#define BASE double
#include "source_tpsv_r.h"
#undef BASE
}

void
cblas_ctpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *Ap, void *X, const int incX)
{
#define BASE float
#include "source_tpsv_c.h"
#undef BASE
}

void
cblas_ztpsv (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
	     const int N, const void *Ap, void *X, const int incX)
{
#define BASE double
#include "source_tpsv_c.h"
#undef BASE
}



/* GER */

void
cblas_sger (const enum CBLAS_ORDER order, const int M, const int N,
	    const float alpha, const float *X, const int incX, const float *Y,
	    const int incY, float *A, const int lda)
{
#define BASE float
#include "source_ger.h"
#undef BASE
}

void
cblas_dger (const enum CBLAS_ORDER order, const int M, const int N,
	    const double alpha, const double *X, const int incX,
	    const double *Y, const int incY, double *A, const int lda)
{
#define BASE double
#include "source_ger.h"
#undef BASE
}


/* GERU */

void
cblas_cgeru (const enum CBLAS_ORDER order, const int M, const int N,
	     const void *alpha, const void *X, const int incX, const void *Y,
	     const int incY, void *A, const int lda)
{
#define BASE float
#include "source_geru.h"
#undef BASE
}

void
cblas_zgeru (const enum CBLAS_ORDER order, const int M, const int N,
	     const void *alpha, const void *X, const int incX, const void *Y,
	     const int incY, void *A, const int lda)
{
#define BASE double
#include "source_geru.h"
#undef BASE
}


/* GERC */

void
cblas_cgerc (const enum CBLAS_ORDER order, const int M, const int N,
	     const void *alpha, const void *X, const int incX, const void *Y,
	     const int incY, void *A, const int lda)
{
#define BASE float
#include "source_gerc.h"
#undef BASE
}

void
cblas_zgerc (const enum CBLAS_ORDER order, const int M, const int N,
	     const void *alpha, const void *X, const int incX, const void *Y,
	     const int incY, void *A, const int lda)
{
#define BASE double
#include "source_gerc.h"
#undef BASE
}

/* HER */

void
cblas_cher (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const float alpha, const void *X, const int incX,
	    void *A, const int lda)
{
#define BASE float
#include "source_her.h"
#undef BASE
}

void
cblas_zher (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const double alpha, const void *X, const int incX,
	    void *A, const int lda)
{
#define BASE double
#include "source_her.h"
#undef BASE
}


/* HPR */

void
cblas_chpr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const float alpha, const void *X, const int incX,
	    void *Ap)
{
#define BASE float
#include "source_hpr.h"
#undef BASE
}

void
cblas_zhpr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const double alpha, const void *X, const int incX,
	    void *Ap)
{
#define BASE double
#include "source_hpr.h"
#undef BASE
}


/* HER2 */

void
cblas_cher2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *X, const int incX,
	     const void *Y, const int incY, void *A, const int lda)
{
#define BASE float
#include "source_her2.h"
#undef BASE
}

void
cblas_zher2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *X, const int incX,
	     const void *Y, const int incY, void *A, const int lda)
{
#define BASE double
#include "source_her2.h"
#undef BASE
}


/* HPR2 */

void
cblas_chpr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *X, const int incX,
	     const void *Y, const int incY, void *Ap)
{
#define BASE float
#include "source_hpr2.h"
#undef BASE
}

void
cblas_zhpr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const void *alpha, const void *X, const int incX,
	     const void *Y, const int incY, void *Ap)
{
#define BASE double
#include "source_hpr2.h"
#undef BASE
}

/* SYR */

void
cblas_ssyr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const float alpha, const float *X, const int incX,
	    float *A, const int lda)
{
#define BASE float
#include "source_syr.h"
#undef BASE
}

void
cblas_dsyr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const double alpha, const double *X, const int incX,
	    double *A, const int lda)
{
#define BASE double
#include "source_syr.h"
#undef BASE
}


/* SPR */

void
cblas_sspr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const float alpha, const float *X, const int incX,
	    float *Ap)
{
#define BASE float
#include "source_spr.h"
#undef BASE
}

void
cblas_dspr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	    const int N, const double alpha, const double *X, const int incX,
	    double *Ap)
{
#define BASE double
#include "source_spr.h"
#undef BASE
}


/* SYR2 */

void
cblas_ssyr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const float alpha, const float *X, const int incX,
	     const float *Y, const int incY, float *A, const int lda)
{
#define BASE float
#include "source_syr2.h"
#undef BASE
}

void
cblas_dsyr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const double alpha, const double *X, const int incX,
	     const double *Y, const int incY, double *A, const int lda)
{
#define BASE double
#include "source_syr2.h"
#undef BASE
}


/* SPR2 */

void
cblas_sspr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const float alpha, const float *X, const int incX,
	     const float *Y, const int incY, float *A)
{
#define BASE double
#include "source_spr2.h"
#undef BASE
}

void
cblas_dspr2 (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
	     const int N, const double alpha, const double *X, const int incX,
	     const double *Y, const int incY, double *A)
{
#define BASE double
#include "source_spr2.h"
#undef BASE
}

/*
 * ===========================================================================
 * level 3 BLAS
 * ===========================================================================
 */

/* GEMM */

void
cblas_sgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	     const int K, const float alpha, const float *A, const int lda,
	     const float *B, const int ldb, const float beta, float *C,
	     const int ldc)
{
#define BASE float
#include "source_gemm_r.h"
#undef BASE
}

void
cblas_dgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	     const int K, const double alpha, const double *A, const int lda,
	     const double *B, const int ldb, const double beta, double *C,
	     const int ldc)
{
#define BASE double
#include "source_gemm_r.h"
#undef BASE
}

void
cblas_cgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	     const int K, const void *alpha, const void *A, const int lda,
	     const void *B, const int ldb, const void *beta, void *C,
	     const int ldc)
{
#define BASE float
#include "source_gemm_c.h"
#undef BASE
}

void
cblas_zgemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	     const int K, const void *alpha, const void *A, const int lda,
	     const void *B, const int ldb, const void *beta, void *C,
	     const int ldc)
{
#define BASE double
#include "source_gemm_c.h"
#undef BASE
}


/* SYMM */

void
cblas_ssymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const int M, const int N,
	     const float alpha, const float *A, const int lda, const float *B,
	     const int ldb, const float beta, float *C, const int ldc)
{
#define BASE float
#include "source_symm_r.h"
#undef BASE
}

void
cblas_dsymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const int M, const int N,
	     const double alpha, const double *A, const int lda,
	     const double *B, const int ldb, const double beta, double *C,
	     const int ldc)
{
#define BASE double
#include "source_symm_r.h"
#undef BASE
}

void
cblas_csymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const int M, const int N,
	     const void *alpha, const void *A, const int lda, const void *B,
	     const int ldb, const void *beta, void *C, const int ldc)
{
#define BASE float
#include "source_symm_c.h"
#undef BASE
}

void
cblas_zsymm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const int M, const int N,
	     const void *alpha, const void *A, const int lda, const void *B,
	     const int ldb, const void *beta, void *C, const int ldc)
{
#define BASE double
#include "source_symm_c.h"
#undef BASE
}


/* HEMM */

void
cblas_chemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	     const void *alpha, const void *A, const int lda, const void *B,
	     const int ldb, const void *beta, void *C, const int ldc)
{
#define BASE float
#include "source_hemm.h"
#undef BASE
}

void
cblas_zhemm (const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
	     const void *alpha, const void *A, const int lda, const void *B,
	     const int ldb, const void *beta, void *C, const int ldc)
{
#define BASE double
#include "source_hemm.h"
#undef BASE
}


/* SYRK */

void
cblas_ssyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	     const float alpha, const float *A, const int lda, float beta,
	     float *C, const int ldc)
{
#define BASE float
#include "source_syrk_r.h"
#undef BASE
}

void
cblas_dsyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	     const double alpha, const double *A, const int lda,
	     const double beta, double *C, const int ldc)
{
#define BASE double
#include "source_syrk_r.h"
#undef BASE
}

void
cblas_csyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	     const void *alpha, const void *A, const int lda,
	     const void *beta, void *C, const int ldc)
{
#define BASE float
#include "source_sryk_c.h"
#undef BASE
}

void
cblas_zsyrk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	     const void *alpha, const void *A, const int lda,
	     const void *beta, void *C, const int ldc)
{
#define BASE double
#include "source_sryk_c.h"
#undef BASE
}


/* HERK */

void
cblas_cherk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	     const float alpha, const void *A, const int lda,
	     const float beta, void *C, const int ldc)
{
#define BASE float
#include "source_herk.h"
#undef BASE
}

void
cblas_zherk (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	     const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	     const double alpha, const void *A, const int lda,
	     const double beta, void *C, const int ldc)
{
#define BASE double
#include "source_herk.h"
#undef BASE
}


/* SYR2K */

void
cblas_ssyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	      const float alpha, const float *A, const int lda,
	      const float *B, const int ldb, const float beta, float *C,
	      const int ldc)
{
#define BASE float
#include "source_syr2k_r.h"
#undef BASE
}

void
cblas_dsyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	      const double alpha, const double *A, const int lda,
	      const double *B, const int ldb, const double beta, double *C,
	      const int ldc)
{
#define BASE double
#include "source_syr2k_r.h"
#undef BASE
}

void
cblas_csyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	      const void *alpha, const void *A, const int lda, const void *B,
	      const int ldb, const void *beta, void *C, const int ldc)
{
#define BASE float
#include "source_syr2k_c.h"
#undef BASE
}

void
cblas_zsyr2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	      const void *alpha, const void *A, const int lda, const void *B,
	      const int ldb, const void *beta, void *C, const int ldc)
{
#define BASE double
#include "source_syr2k_c.h"
#undef BASE
}


/* HER2K */

void
cblas_cher2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	      const void *alpha, const void *A, const int lda, const void *B,
	      const int ldb, const float beta, void *C, const int ldc)
{
#define BASE float
#include "source_her2k.h"
#undef BASE
}


void
cblas_zher2k (const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
	      const enum CBLAS_TRANSPOSE Trans, const int N, const int K,
	      const void *alpha, const void *A, const int lda, const void *B,
	      const int ldb, const double beta, void *C, const int ldc)
{
#define BASE double
#include "source_her2k.h"
#undef BASE
}


/* TRMM */

void
cblas_strmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const float alpha, const float *A, const int lda, float *B,
	     const int ldb)
{
#define BASE float
#include "source_trmm_r.h"
#undef BASE
}

void
cblas_dtrmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const double alpha, const double *A, const int lda, double *B,
	     const int ldb)
{
#define BASE double
#include "source_trmm_r.h"
#undef BASE

}

void
cblas_ctrmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const void *alpha, const void *A, const int lda, void *B,
	     const int ldb)
{
#define BASE float
#include "source_trmm_c.h"
#undef BASE
}

void
cblas_ztrmm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const void *alpha, const void *A, const int lda, void *B,
	     const int ldb)
{
#define BASE double
#include "source_trmm_c.h"
#undef BASE
}


/* TRSM */

void
cblas_strsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const float alpha, const float *A, const int lda, float *B,
	     const int ldb)
{
#define BASE float
#include "source_trsm_r.h"
#undef BASE

}

void
cblas_dtrsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const double alpha, const double *A, const int lda, double *B,
	     const int ldb)
{
#define BASE double
#include "source_trsm_r.h"
#undef BASE
}

void
cblas_ctrsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const void *alpha, const void *A, const int lda, void *B,
	     const int ldb)
{
#define BASE float
#include "source_trsm_c.h"
#undef BASE
}

void
cblas_ztrsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
	     const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
	     const enum CBLAS_DIAG Diag, const int M, const int N,
	     const void *alpha, const void *A, const int lda, void *B,
	     const int ldb)
{
#define BASE double
#include "source_trsm_c.h"
#undef BASE
}
#endif
