/*
 * Author:  G. Jungman  and  B. Gough
 * RCS:     $Id$
 */
/* Native implementation of row-major operations.
 * Conforms to gsl_blas_raw interface.
 */
#include <math.h>
#include "gsl_blas_raw.h"


/* we need some local macros for handling the "raw" complex data */

/* void * X :=  float * X */
#define S_REAL(X, i)  (*((float *)X + 2*i))
#define S_IMAG(X, i)  (*((float *)X + 2*i + 1))

/* void * X :=  double * X */
#define D_REAL(X, i)  (*((double *)X + 2*i))
#define D_IMAG(X, i)  (*((double *)X + 2*i + 1))



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
                         const void * X, int incX,
                         const void * Y, int incY,
                         void * dotu)
{
  float rr = 0.0;
  float ri = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += S_REAL(X, i)*S_REAL(Y, j) - S_IMAG(X, i)*S_IMAG(Y, j);
    ri += S_REAL(X, i)*S_IMAG(Y, j) + S_IMAG(X, i)*S_REAL(Y, j);
    i += incX;
    j += incY;
  }
  S_REAL(dotu, 0) = rr;
  S_REAL(dotu, 0) = ri;
}

void gsl_blas_raw_cdotc (size_t N,
                         const void * X, int incX,
                         const void * Y, int incY,
                         void * dotc)
{
  float rr = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += S_REAL(X, i)*S_REAL(Y, j) + S_IMAG(X, i)*S_IMAG(Y, j);
    i += incX;
    j += incY;
  }
  S_REAL(dotc, 0) = rr;
  S_IMAG(dotc, 0) = 0.0;
}

void gsl_blas_raw_zdotu (size_t N,
                         const void * X, int incX,
                         const void * Y, int incY,
                         void * dotu)
{
  float rr = 0.0;
  float ri = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += D_REAL(X, i)*D_REAL(Y, j) - D_IMAG(X, i)*D_IMAG(Y, j);
    ri += D_REAL(X, i)*D_IMAG(Y, j) + D_IMAG(X, i)*D_REAL(Y, j);
    i += incX;
    j += incY;
  }
  D_REAL(dotu, 0) = rr;
  D_IMAG(dotu, 0) = ri; 
}

void gsl_blas_raw_zdotc (size_t N,
                         const void * X, int incX,
                         const void * Y, int incY,
                         void * dotc)
{
  double rr = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += D_REAL(X, i)*D_REAL(Y, j) + D_IMAG(X, i)*D_IMAG(Y, j);
    i += incX;
    j += incY;
  }
  D_REAL(dotc, 0) = rr;
  D_IMAG(dotc, 0) = 0.0;
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

float  gsl_blas_raw_scnrm2 (size_t N, const void * X, int incX)
{
  float r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += S_REAL(X, i)*S_REAL(X, i) + S_IMAG(X, i)*S_IMAG(X, i);
    i += incX;
  }
  return sqrt(r);
}

double gsl_blas_raw_dznrm2 (size_t N, const void * X, int incX)
{
  double r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += D_REAL(X, i)*D_REAL(X, i) + D_IMAG(X, i)*D_IMAG(X, i);
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

float  gsl_blas_raw_scasum (size_t N, const void * X, int incX)
{
  float r = 0.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += fabs(S_REAL(X, i)) + fabs(S_IMAG(X, i));
    i += incX;
  }
  return r;
}

double gsl_blas_raw_dzasum (size_t N, const void * X, int incX)
{
 double r = 0.0;
 size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    r += fabs(D_REAL(X, i)) + fabs(D_IMAG(X, i));
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

CBLAS_INDEX gsl_blas_raw_icamax (size_t N, const void * X, int incX)
{
  double max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    double a = fabs(S_REAL(X, i)) + fabs(S_IMAG(X, i));
    if(a > max) {
      max = a;
      result = i;
    }
    i += incX;
  }
  return result;
}

CBLAS_INDEX gsl_blas_raw_izamax (size_t N, const void * X, int incX)
{
  double max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    double a = fabs(D_REAL(X, i)) + fabs(D_IMAG(X, i));
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
                         void * X, int incX,
                         void * Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    float tmpr = S_REAL(X, i);
    float tmpi = S_IMAG(X, i);
    S_REAL(X, i) = S_REAL(Y, j);
    S_IMAG(X, i) = S_IMAG(Y, j);
    S_REAL(Y, j) = tmpr;
    S_IMAG(Y, j) = tmpi;
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_zswap (size_t N,
                         void * X, int incX,
                         void * Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    float tmpr = D_REAL(X, i);
    float tmpi = D_IMAG(X, i);
    D_REAL(X, i) = D_REAL(Y, j);
    D_IMAG(X, i) = D_IMAG(Y, j);
    D_REAL(Y, j) = tmpr;
    D_IMAG(Y, j) = tmpi;
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
                            const void * X, int incX,
                            void * Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    S_REAL(Y, j) = S_REAL(X, i);
    S_IMAG(Y, j) = S_IMAG(X, i);
    i += incX;
    j += incY;
  }
}

void gsl_blas_raw_zcopy (size_t N,
                         const void * X, int incX,
                         void * Y, int incY)
{
  size_t n;
  size_t i;
  size_t j;
  X = (gsl_complex *) X;
  Y = (gsl_complex *) Y;
  for(n=0; n<N; n++) {
    D_REAL(Y, j) = D_REAL(X, i);
    D_IMAG(Y, j) = D_IMAG(X, i);
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
                         const void * alpha,
                         const void * X, int incX,
                         void * Y, int incY)
{
}

void gsl_blas_raw_zaxpy (size_t N,
                         const void * alpha,
                         const void * X, int incX,
                         void * Y, int incY)
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
void gsl_blas_raw_cscal  (size_t N, const void * alpha, void * X, int incX)
{
}
void gsl_blas_raw_zscal  (size_t N, const void * alpha, void * X, int incX)
{
}
void gsl_blas_raw_csscal (size_t N, float  alpha, void * X, int incX)
{
}
void gsl_blas_raw_zdscal (size_t N, double alpha, void * X, int incX)
{
}


/* ===========================================================================
 * Generate the row-major versions of the level 2 and level 3 operations.
 * ===========================================================================
 */

#define  MACCESS(s, i, j)  ((s)*(i) + (j))
#define  FUNC(f) f
#include "blas_raw_native_L23_source.c"
#undef   FUNC
#undef   MACCESS

