/*
 * Author:  G. Jungman  and  B. Gough
 * RCS:     $Id$
 */
/* Native implementation of row-major operations.
 * Conforms to gsl_blas_raw interface.
 */
#include <math.h>
#include "complex_internal.h"
#include "gsl_blas_raw.h"



/* ===========================================================================
 * level 1 BLAS functions
 * ===========================================================================
 */

float  gsl_blas_raw_sdsdot (size_t N,
                            float alpha,
                            const float X[], int incX,
                            const float Y[], int incY)
{ /* FIXME: is this right ?? */
  float r = 0.0;
  size_t n;
  size_t i = 0;
  size_t j = 0;
  for(n=0; n<N; n++) {
    r += X[i]*Y[j];
    i += incX;
    j += incY;
  }
  return r + alpha;
}

double gsl_blas_raw_dsdot (size_t N,
                           const float X[], int incX,
                           const float Y[], int incY)
{
  double r = 0.0;
  size_t n;
  size_t i = 0;
  size_t j = 0;
  for(n=0; n<N; n++) {
    r += X[i]*Y[j];
    i += incX;
    j += incY;
  }
  return r;
}

float  gsl_blas_raw_sdot (size_t N,
                          const float X[], int incX,
                          const float Y[], int incY)
{
  float r = 0.0;
  size_t n;
  size_t i = 0;
  size_t j = 0;
  for(n=0; n<N; n++) {
    r += X[i]*Y[j];
    i += incX;
    j += incY;
  }
  return r;
}
                          
double gsl_blas_raw_ddot (size_t N,
                          const double X[], int incX,
                          const double Y[], int incY)
{
  double r = 0.0;
  size_t n;
  size_t i = 0;
  size_t j = 0;
  for(n=0; n<N; n++) {
    r += X[i]*Y[j];
    i += incX;
    j += incY;
  }
  return r;
}


void gsl_blas_raw_cdotu (size_t N,
                         const gsl_complex_packed_array_float X, int incX,
                         const gsl_complex_packed_array_float Y, int incY,
                         gsl_complex_packed_float dotu)
{
  float rr = 0.0;
  float ri = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += REAL(X, incX, i)*REAL(Y, incY, j) - IMAG(X, incX, i)*IMAG(Y, incY, j);
    ri += REAL(X, incX, i)*IMAG(Y, incY, j) + IMAG(X, incX, i)*REAL(Y, incY, j);
    i += incX;
    j += incY;
  }
  REAL0(dotu) = rr;
  IMAG0(dotu) = ri;
}

void gsl_blas_raw_cdotc (size_t N,
                         const gsl_complex_packed_array_float X, int incX,
                         const gsl_complex_packed_array_float Y, int incY,
                         gsl_complex_packed_float dotc)
{
  float rr = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += REAL(X, incX, i)*REAL(Y, incY, j) + IMAG(X, incX, i)*IMAG(Y, incY, j);
    i += incX;
    j += incY;
  }
  REAL0(dotc) = rr;
  IMAG0(dotc) = 0.0;
}

void gsl_blas_raw_zdotu (size_t N,
                         const gsl_complex_packed_array X, int incX,
                         const gsl_complex_packed_array Y, int incY,
                         gsl_complex_packed dotu)
{
  float rr = 0.0;
  float ri = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += REAL(X, incX, i)*REAL(Y, incY, j) - IMAG(X, incX, i)*IMAG(Y, incY, j);
    ri += REAL(X, incX, i)*IMAG(Y, incY, j) + IMAG(X, incX, i)*REAL(Y, incY, j);
    i += incX;
    j += incY;
  }
  REAL0(dotu) = rr;
  IMAG0(dotu) = ri; 
}

void gsl_blas_raw_zdotc (size_t N,
                         const gsl_complex_packed_array X, int incX,
                         const gsl_complex_packed_array Y, int incY,
                         gsl_complex_packed dotc)
{
  double rr = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += REAL(X, incX, i)*REAL(Y, incY, j) + IMAG(X, incX, i)*IMAG(Y, incY, j);
    i += incX;
    j += incY;
  }
  REAL0(dotc) = rr;
  IMAG0(dotc) = 0.0;
}


float  gsl_blas_raw_snrm2  (size_t N, const float  X[], int incX)
{
  float r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += X[i]*X[i];
    i += incX;
  }
  return sqrt(r);
}

double gsl_blas_raw_dnrm2  (size_t N, const double X[], int incX)
{
  double r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += X[i]*X[i];
    i += incX;
  }
  return sqrt(r);
}

float  gsl_blas_raw_scnrm2 (size_t N, const gsl_complex_packed_array_float X, int incX)
{
  float r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += REAL(X, incX, i)*REAL(X, incX, i) + IMAG(X, incX, i)*IMAG(X, incX, i);
    i += incX;
  }
  return sqrt(r);
}

double gsl_blas_raw_dznrm2 (size_t N, const gsl_complex_packed_array X, int incX)
{
  double r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += REAL(X, incX, i)*REAL(X, incX, i) + IMAG(X, incX, i)*IMAG(X, incX, i);
    i += incX;
  }
  return sqrt(r);
}

float  gsl_blas_raw_sasum (size_t N, const float X[], int incX)
{
  float r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += fabs(X[i]);
    i += incX;
  }
  return r;
}

double gsl_blas_raw_dasum (size_t N, const double X[], int incX)
{
  double r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += fabs(X[i]);
    i += incX;
  }
  return r;
}

float  gsl_blas_raw_scasum (size_t N, const gsl_complex_packed_array_float X, int incX)
{
  float r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    i += incX;
  }
  return r;
}

double gsl_blas_raw_dzasum (size_t N, const gsl_complex_packed_array X, int incX)
{
 double r = 0.0;
 size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    i += incX;
  }
 return r;
}


CBLAS_INDEX gsl_blas_raw_isamax (size_t N, const float X[], int incX)
{
  float max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    if(fabs(X[i]) > max) {
      max = fabs(X[i]);
      result = i;
    }
    i += incX;
  }
  return result;
}

CBLAS_INDEX gsl_blas_raw_idamax (size_t N, const double X[], int incX)
{
  double max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    if(fabs(X[i]) > max) {
      max = fabs(X[i]);
      result = i;
    }
    i += incX;
  }
  return result;
}

CBLAS_INDEX gsl_blas_raw_icamax (size_t N, const gsl_complex_packed_array_float X, int incX)
{
  double max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    double a = fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    if(a > max) {
      max = a;
      result = i;
    }
    i += incX;
  }
  return result;
}

CBLAS_INDEX gsl_blas_raw_izamax (size_t N, const gsl_complex_packed_array X, int incX)
{
  double max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    double a = fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    if(a > max) {
      max = a;
      result = i;
    }
    i += incX;
  }
  return result;
}


void gsl_blas_raw_sswap (size_t N,
                         float X[], int incX,
                         float Y[], int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    double tmp = X[i];
    X[i] = Y[j];
    Y[j] = tmp;
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_dswap (size_t N,
                         double X[], int incX,
                         double Y[], int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    double tmp = X[i];
    X[i] = Y[j];
    Y[j] = tmp;
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_cswap (size_t N,
                         gsl_complex_packed_array_float X, int incX,
                         gsl_complex_packed_array_float Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    float tmpr = REAL(X, incX, i);
    float tmpi = IMAG(X, incX, i);
    REAL(X, incX, i) = REAL(Y, incY, j);
    IMAG(X, incX, i) = IMAG(Y, incY, j);
    REAL(Y, incY, j) = tmpr;
    IMAG(Y, incY, j) = tmpi;
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_zswap (size_t N,
                         gsl_complex_packed_array X, int incX,
                         gsl_complex_packed_array Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    float tmpr = REAL(X, incX, i);
    float tmpi = IMAG(X, incX, i);
    REAL(X, incX, i) = REAL(Y, incY, j);
    IMAG(X, incX, i) = IMAG(Y, incY, j);
    REAL(Y, incY, j) = tmpr;
    IMAG(Y, incY, j) = tmpi;
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_scopy (size_t N,
                         const float X[], int incX,
                         float Y[], int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    Y[j] = X[i];
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_dcopy (size_t N,
                         const double X[], int incX,
                         double Y[], int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    Y[j] = X[i];
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_ccopy (size_t N,
                         const gsl_complex_packed_array_float X, int incX,
                         gsl_complex_packed_array_float Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    REAL(Y, incY, j) = REAL(X, incX, i);
    IMAG(Y, incY, j) = IMAG(X, incX, i);
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_zcopy (size_t N,
                         const gsl_complex_packed_array X, int incX,
                         gsl_complex_packed_array Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    REAL(Y, incY, j) = REAL(X, incX, i);
    IMAG(Y, incY, j) = IMAG(X, incX, i);
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_saxpy (size_t N,
                         float alpha,
                         const float X[], int incX,
                         float Y[], int incY)
{
}


void gsl_blas_raw_daxpy (size_t N,
                         double alpha,
                         const double X[], int incX, 
                         double Y[], int incY)
{
}


void gsl_blas_raw_caxpy (size_t N,
                         const gsl_complex_packed_float alpha,
                         const gsl_complex_packed_array_float X, int incX,
                         gsl_complex_packed_array_float Y, int incY)
{
}

void gsl_blas_raw_zaxpy (size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array X, int incX,
                         gsl_complex_packed_array Y, int incY)
{
}


void gsl_blas_raw_srotg (float a[], float b[], float c[], float s[])
{
}
void gsl_blas_raw_drotg (double a[], double b[], double c[], double s[])
{
}

void gsl_blas_raw_srotmg (float d1[], float d2[], float b1[], float b2, float P[])
{
}
void gsl_blas_raw_drotmg (double d1[], double d2[], double b1[],
                             double b2, double P[])
{
}


void gsl_blas_raw_srot (size_t N,
                           float X[], int incX,
                           float Y[], int incY,
                           float c, float s)
{
}

void gsl_blas_raw_drot (size_t N,
                           double X[], int incX,
                           double Y[], int incY,
                           const double c, const double s)
{
}

void gsl_blas_raw_srotm (size_t N,
                            float X[], int incX,
                            float Y[], int incY,
                            const float P[])
{
}

void gsl_blas_raw_drotm (size_t N,
                            double X[], int incX,
                            double Y[], int incY,
                            const double P[])
{
}


void gsl_blas_raw_sscal  (size_t N, float  alpha, float  X[], int incX)
{
}
void gsl_blas_raw_dscal  (size_t N, double alpha, double X[], int incX)
{
}
void gsl_blas_raw_cscal  (size_t N, const gsl_complex_packed_float alpha, gsl_complex_packed_array_float X, int incX)
{
}
void gsl_blas_raw_zscal  (size_t N, const gsl_complex_packed alpha, gsl_complex_packed_array X, int incX)
{
}
void gsl_blas_raw_csscal (size_t N, float  alpha, gsl_complex_packed_array_float X, int incX)
{
}
void gsl_blas_raw_zdscal (size_t N, double alpha, gsl_complex_packed_array X, int incX)
{
}


/*
 * ===========================================================================
 * level 2 BLAS
 * ===========================================================================
 */
 
#define FUNC(f) f


/* GEMV */

void FUNC(gsl_blas_raw_sgemv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N,
                               float alpha,
                               const float A[], int lda,
                               const float X[], int incX,
                               float beta,
                               float Y[], int incY)
{
  size_t i, j;

  if(TransA == CblasNoTrans) {
    for(i=0; i<M; i++) {
      for(j=0; j<N; j++) {
        Y[i * incY] = alpha * A[MACCESS(lda, i ,j)] * X[j * incX]
	              + beta * Y[i * incY];
      }
    }
  }
  else {  /* CBlasTrans || CBlasConjTrans */
    for(i=0; i<N; i++) {
      for(j=0; j<M; j++) {
        Y[i * incY] = alpha * A[MACCESS(lda, j ,i)] * X[j * incX]
	              + beta * Y[i * incY];
      }
    }
  }
}

void FUNC(gsl_blas_raw_dgemv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N,
                               double alpha,
                               const double A[], int lda,
                               const double X[], int incX,
                               double beta,
                               double Y[], int incY)
{
}

void FUNC(gsl_blas_raw_cgemv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float Y, int incY)
{
}

void FUNC(gsl_blas_raw_zgemv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, int incY)
{
}

/* GBMV */

void FUNC(gsl_blas_raw_sgbmv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N, size_t KL, size_t KU,
                               float alpha,
                               const float A[], int lda,
                               const float X[], int incX,
                               float beta,
                               float Y[], int incY)
{
}

void FUNC(gsl_blas_raw_dgbmv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N, size_t KL, size_t KU,
                               double alpha,
                               const double A[], int lda,
                               const double X[], int incX,
                               double beta,
                               double Y[], int incY)
{
}

void FUNC(gsl_blas_raw_cgbmv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N, size_t KL, size_t KU,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float Y, int incY)
{
}

void FUNC(gsl_blas_raw_zgbmv) (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N, size_t KL, size_t KU,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, int incY)
{
}


/* TRMV */

void FUNC(gsl_blas_raw_strmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const float A[], int lda,
                               float X[], int incX)
{
}

void FUNC(gsl_blas_raw_dtrmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double A[], int lda,
                               double X[], int incX)
{
}

void FUNC(gsl_blas_raw_ctrmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, int incX)
{
}

void FUNC(gsl_blas_raw_ztrmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, int incX)
{
}


/* TBMV */

void FUNC(gsl_blas_raw_stbmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const float A[], int lda,
                               float X[], int incX)
{
}

void FUNC(gsl_blas_raw_dtbmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const double A[], int lda,
                               double X[], int incX)
{
}

void FUNC(gsl_blas_raw_ctbmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, int incX)
{
}

void FUNC(gsl_blas_raw_ztbmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, int incX)
{
}


/* TPMV */

void FUNC(gsl_blas_raw_stpmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const float Ap[],
                               float X[], int incX)
{
}

void FUNC(gsl_blas_raw_dtpmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double Ap[],
                               double X[], int incX)
{
}

void FUNC(gsl_blas_raw_ctpmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float Ap,
                               gsl_complex_packed_array_float X, int incX)
{
}

void FUNC(gsl_blas_raw_ztpmv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array Ap,
                               gsl_complex_packed_array X, int incX)
{
}

/* TRSV */

void FUNC(gsl_blas_raw_strsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const float A[], int lda,
                               float X[], int incX)
{
}

void FUNC(gsl_blas_raw_dtrsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double A[], int lda,
                               double X[], int incX)
{
}

void FUNC(gsl_blas_raw_ctrsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, int incX)
{
}

void FUNC(gsl_blas_raw_ztrsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, int incX)
{
}


/* TBSV */

void FUNC(gsl_blas_raw_stbsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const float A[], int lda,
                               float X[], int incX)
{
}

void FUNC(gsl_blas_raw_dtbsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const double A[], int lda,
                               double X[], int incX)
{
}

void FUNC(gsl_blas_raw_ctbsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, int incX)
{
}

void FUNC(gsl_blas_raw_ztbsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, int incX)
{
}


/* TPSV */

void FUNC(gsl_blas_raw_stpsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const float Ap[],
                               float X[], int incX)
{
}

void FUNC(gsl_blas_raw_dtpsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double Ap[],
                               double X[], int incX)
{
}

void FUNC(gsl_blas_raw_ctpsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float Ap,
                               gsl_complex_packed_array_float X, int incX)
{
}

void FUNC(gsl_blas_raw_ztpsv) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array Ap,
                               gsl_complex_packed_array X, int incX)
{
}


/* SYMV */

void FUNC(gsl_blas_raw_ssymv) (CBLAS_UPLO Uplo,
                               size_t N,
                               float alpha,
                               const float A[], int lda,
                               const float X[], int incX,
                               float beta,
                               float Y[], int incY)
{
}

void FUNC(gsl_blas_raw_dsymv) (CBLAS_UPLO Uplo,
                               size_t N,
                               double alpha,
                               const double A[], int lda,
                               const double X[], int incX,
                               double beta,
                               double Y[], int incY)
{
}


/* SBMV */

void FUNC(gsl_blas_raw_ssbmv) (CBLAS_UPLO Uplo,
                               size_t N, size_t K,
                               float alpha,
                               const float A[], int lda, 
                               const float X[], int incX,
                               float beta,
                               float Y[], int incY)
{
}

void FUNC(gsl_blas_raw_dsbmv) (CBLAS_UPLO Uplo,
                               size_t N, size_t K,
                               double alpha,
                               const double A[], int lda,
                               const double X[], int incX,
                               double beta,
                               double Y[], int incY)
{
}

/* SPMV */

void FUNC(gsl_blas_raw_sspmv) (CBLAS_UPLO Uplo,
                               size_t N,
                               float alpha,
                               const float Ap[],
                               const float X[], int incX,
                               float beta,
                               float Y[], int incY)
{
}

void FUNC(gsl_blas_raw_dspmv) (CBLAS_UPLO Uplo,
                               size_t N,
                               double alpha,
                               const double Ap[],
                               const double X[], int incX,
                               double beta,
                               double Y[], int incY)
{
}

/* GER */

void FUNC(gsl_blas_raw_sger) (size_t M, size_t N,
                              float alpha,
                              const float X[], int incX,
                              const float Y[], int incY,
                              float A[], int lda)
{
}

void FUNC(gsl_blas_raw_dger) (size_t M, size_t N,
                              double alpha,
                              const double X[], int incX,
                              const double Y[], int incY,
                              double A[], int lda)
{
}


/* SYR */

void FUNC(gsl_blas_raw_ssyr) (CBLAS_UPLO Uplo,
                              size_t N,
                              float alpha,
                              const float X[], int incX,
                              float A[], int lda)
{
}

void FUNC(gsl_blas_raw_dsyr) (CBLAS_UPLO Uplo,
                              size_t N,
                              double alpha,
                              const double X[], int incX,
                              double A[], int lda)
{
}


/* SPR */

void FUNC(gsl_blas_raw_sspr) (CBLAS_UPLO Uplo,
                              size_t N,
                              float alpha,
                              const float X[], int incX,
                              float Ap[])
{
}

void FUNC(gsl_blas_raw_dspr) (CBLAS_UPLO Uplo,
                              size_t N,
                              double alpha,
                              const double X[], int incX,
                              double Ap[])
{
}


/* SYR2 */

void FUNC(gsl_blas_raw_ssyr2) (CBLAS_UPLO Uplo,
                               size_t N,
                               float alpha,
                               const float X[], int incX, 
                               const float Y[], int incY,
                               float A[], int lda)
{
}

void FUNC(gsl_blas_raw_dsyr2) (CBLAS_UPLO Uplo,
                               size_t N,
                               double alpha,
                               const double X[], int incX,
                               const double Y[], int incY,
                               double A[], int lda)
{
}


/* SPR2 */

void FUNC(gsl_blas_raw_sspr2) (CBLAS_UPLO Uplo,
                               size_t N,
                               float alpha,
                               const float X[], int incX,
                               const float Y[], int incY,
                               float A[])
{
}

void FUNC(gsl_blas_raw_dspr2) (CBLAS_UPLO Uplo,
                               size_t N,
                               double alpha,
                               const double X[], int incX,
                               const double Y[], int incY,
                               double A[])
{
}


/* HEMV */

void FUNC(gsl_blas_raw_chemv) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float Y, int incY)
{
}

void FUNC(gsl_blas_raw_zhemv) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, int incY)
{
}


/* HBMV */

void FUNC(gsl_blas_raw_chbmv) (CBLAS_UPLO Uplo,
                               size_t N, size_t K,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float Y, int incY)
{
}

void FUNC(gsl_blas_raw_zhbmv) (CBLAS_UPLO Uplo,
                               size_t N, size_t K,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, int incY)
{
}


/* HPMV */

void FUNC(gsl_blas_raw_chpmv) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float Ap,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float Y, int incY)
{
}

void FUNC(gsl_blas_raw_zhpmv) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array Ap,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, int incY)
{
}


/* GERU */

void FUNC(gsl_blas_raw_cgeru) (size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_array_float Y, int incY,
                               gsl_complex_packed_array_float A, int lda)
{
}

void FUNC(gsl_blas_raw_zgeru) (size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed_array Y, int incY,
                               gsl_complex_packed_array A, int lda)
{
}


/* GERC */

void FUNC(gsl_blas_raw_cgerc) (size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_array_float Y, int incY,
                               gsl_complex_packed_array_float A, int lda)
{
}

void FUNC(gsl_blas_raw_zgerc) (size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed_array Y, int incY,
                               gsl_complex_packed_array A, int lda)
{
}

/* HER */

void FUNC(gsl_blas_raw_cher) (CBLAS_UPLO Uplo,
                              size_t N,
                              float alpha,
                              const gsl_complex_packed_array_float X, int incX,
                              gsl_complex_packed_array_float A, int lda)
{
}

void FUNC(gsl_blas_raw_zher) (CBLAS_UPLO Uplo,
                              size_t N,
                              double alpha,
                              const gsl_complex_packed_array X, int incX,
                              gsl_complex_packed_array A, int lda)
{
}


/* HPR */

void FUNC(gsl_blas_raw_chpr) (CBLAS_UPLO Uplo,
                              size_t N,
                              float alpha,
                              const gsl_complex_packed_array_float X, int incX,
                              gsl_complex_packed_array_float A)
{
}

void FUNC(gsl_blas_raw_zhpr) (CBLAS_UPLO Uplo,
                              size_t N,
                              double alpha,
                              const gsl_complex_packed_array X, int incX,
                              gsl_complex_packed_array A)
{
}


/* HER2 */

void FUNC(gsl_blas_raw_cher2) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_array_float Y, int incY,
                               gsl_complex_packed_array_float A, int lda)
{
}

void FUNC(gsl_blas_raw_zher2) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed_array Y, int incY,
                               gsl_complex_packed_array A, int lda)
{
}


/* HPR2 */

void FUNC(gsl_blas_raw_chpr2) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, int incX,
                               const gsl_complex_packed_array_float Y, int incY,
                               gsl_complex_packed_array_float Ap)
{
}

void FUNC(gsl_blas_raw_zhpr2) (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, int incX,
                               const gsl_complex_packed_array Y, int incY,
                               gsl_complex_packed_array Ap)
{
}


/*
 * ===========================================================================
 * level 3 BLAS
 * ===========================================================================
 */

/* GEMM */

void FUNC(gsl_blas_raw_sgemm) (CBLAS_TRANSPOSE TransA,
                               CBLAS_TRANSPOSE TransB,
                               size_t M, size_t N, size_t K,
                               float alpha,
                               const float A[], int lda,
                               const float B[], int ldb,
                               float beta,
                               float C[], int ldc)
{
}

void FUNC(gsl_blas_raw_dgemm) (CBLAS_TRANSPOSE TransA,
                               CBLAS_TRANSPOSE TransB,
                               size_t M, size_t N, size_t K,
                               double alpha,
                               const double A[], int lda,
                               const double B[], int ldb,
                               double beta,
                               double C[], int ldc)
{
}

void FUNC(gsl_blas_raw_cgemm) (CBLAS_TRANSPOSE TransA,
                               CBLAS_TRANSPOSE TransB,
                               size_t M, size_t N, size_t K,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda, 
                               const gsl_complex_packed_array_float B, int ldb,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float C, int ldc)
{
}

void FUNC(gsl_blas_raw_zgemm) (CBLAS_TRANSPOSE TransA,
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

void FUNC(gsl_blas_raw_ssymm) (CBLAS_SIDE Side, CBLAS_UPLO Uplo,
                               size_t M, size_t N,
                               float alpha,
                               const float A[], int lda,
                               const float B[], int ldb,
                               float beta,
                               float C[], int ldc)
{
}

void FUNC(gsl_blas_raw_dsymm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo,
                               size_t M, size_t N,
                               double alpha,
                               const double A[], int lda,
                               const double B[], int ldb,
                               double beta,
                               double C[], int ldc)
{
}

void FUNC(gsl_blas_raw_csymm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float B, int ldb,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float C, int ldc)
{
}

void FUNC(gsl_blas_raw_zsymm) (CBLAS_SIDE Side,
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

void FUNC(gsl_blas_raw_ssyrk) (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               float alpha,
                               const float A[], int lda,
                               float beta,
                               float C[], int ldc)
{
}

void FUNC(gsl_blas_raw_dsyrk) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               double alpha,
                               const double A[], int lda,
                               double beta,
                               double C[], int ldc)
{
}

void FUNC(gsl_blas_raw_csyrk) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float C, int ldc)
{
}

void FUNC(gsl_blas_raw_zsyrk) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array C, int ldc)
{
}


/* SYR2K */

void FUNC(gsl_blas_raw_ssyr2k) (CBLAS_UPLO Uplo, CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                float alpha,
                                const float A[], int lda,
                                const float B[], int ldb,
                                float beta,
                                float C[], int ldc)
{
}

void FUNC(gsl_blas_raw_dsyr2k) (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                double alpha,
                                const double A[], int lda,
                                const double B[], int ldb,
                                double beta,
                                double C[], int ldc)
{
}

void FUNC(gsl_blas_raw_csyr2k) (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                const gsl_complex_packed_float alpha,
                                const gsl_complex_packed_array_float A, int lda,
                                const gsl_complex_packed_array_float B, int ldb,
                                const gsl_complex_packed_float beta,
                                gsl_complex_packed_array_float C, int ldc)
{
}

void FUNC(gsl_blas_raw_zsyr2k) (CBLAS_UPLO Uplo,
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

void FUNC(gsl_blas_raw_strmm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               float alpha,
                               const float A[], int lda,
                               float B[], int ldb)
{
}

void FUNC(gsl_blas_raw_dtrmm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               double alpha,
                               const double A[], int lda,
                               double B[], int ldb)
{
}

void FUNC(gsl_blas_raw_ctrmm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float B, int ldb)
{
}

void FUNC(gsl_blas_raw_ztrmm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array B, int ldb)
{
}


/* TRSM */

void FUNC(gsl_blas_raw_strsm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               float alpha,
                               const float A[], int lda,
                               float B[], int ldb)
{
}

void FUNC(gsl_blas_raw_dtrsm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               double alpha,
                               const double A[], int lda,
                               double B[], int ldb)
{
}

void FUNC(gsl_blas_raw_ctrsm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float B, int ldb)
{
}

void FUNC(gsl_blas_raw_ztrsm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array B, int ldb)
{
}


/* HEMM */

void FUNC(gsl_blas_raw_chemm) (CBLAS_SIDE Side,
                               CBLAS_UPLO Uplo,
                               size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float B, int ldb,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float C, int ldc)
{
}

void FUNC(gsl_blas_raw_zhemm) (CBLAS_SIDE Side,
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

void FUNC(gsl_blas_raw_cherk) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               float beta,
                               gsl_complex_packed_array_float C, int ldc)
{
}

void FUNC(gsl_blas_raw_zherk) (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE Trans,
                               size_t N, size_t K,
                               double alpha,
                               const gsl_complex_packed_array A, int lda,
                               double beta,
                               gsl_complex_packed_array C, int ldc)
{
}


/* HER2K */

void FUNC(gsl_blas_raw_cher2k) (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                const gsl_complex_packed_float alpha,
                                const gsl_complex_packed_array_float A, int lda,
                                const gsl_complex_packed_array_float B, int ldb,
                                float beta,
                                gsl_complex_packed_array_float C, int ldc)
{
}


void FUNC(gsl_blas_raw_zher2k) (CBLAS_UPLO Uplo,
                                CBLAS_TRANSPOSE Trans,
                                size_t N, size_t K,
                                const gsl_complex_packed alpha,
                                const gsl_complex_packed_array A, int lda,
                                const gsl_complex_packed_array B, int ldb,
                                double beta,
                                gsl_complex_packed_array C, int ldc)
{
}

