/*
 * Author:  G. Jungman  and  B. Gough
 * RCS:     $Id$
 */
/* Native implementation of row-major operations.
 * Conforms to gsl_blas_raw interface.
 *
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
/*
  float r = alpha;
  size_t n;
  size_t i = 0;
  size_t j = 0;
  for(n=0; n<N; n++) {
    r += X[i]*Y[j];
    i += incX;
    j += incY;
  }
  return r;
  */
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
/*
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
  */
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
/*
  float rr = 0.0;
  float ri = 0.0;
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    rr += REAL(X, incX, i)*REAL(Y, incY, j) + IMAG(X, incX, i)*IMAG(Y, incY, j);
    ri += REAL(X, incX, i)*IMAG(Y, incY, j) - IMAG(X, incX, i)*REAL(Y, incY, j);
    i += incX;
    j += incY;
  }
  REAL0(dotc) = rr;
  IMAG0(dotc) = ri;
  */
}

void gsl_blas_raw_zdotu (size_t N,
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
/*
  float scale = 0.0;
  float ssq   = 1.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    float axi = fabs(X[i]);
    if(scale < axi) {
      ssq   = 1.0 + ssq*(scale/axi)*(scale/axi);
      scale = axi;
    }
    else {
      ssq += (axi/scale)*(axi/scale);
    }
    i += incX;
  }
  return scale * sqrt(ssq);
  */
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
/*
  float scale = 0.0;
  float ssq   = 1.0;
  size_t n;
  size_t i;
  for(n=0; n<N; n++) {
    float axi = fabs(REAL(X,incX,i));
    float ayi = fabs(IMAG(X,incX,i));
    if(scale < axi) {
      ssq   = 1.0 + ssq*(scale/axi)*(scale/axi);
      scale = axi;
    }
    else {
      ssq += (axi/scale)*(axi/scale);
    }
    if(scale < ayi) {
      ssq   = 1.0 + ssq*(scale/ayi)*(scale/ayi);
      scale = axi;
    }
    else {
      ssq += (ayi/scale)*(ayi/scale);
    }
    i += incX;
  }
  return scale * sqrt(ssq);
  */
}

double gsl_blas_raw_dznrm2 (size_t N, const gsl_complex_packed_array X, size_t incX)
{
#define BASE_TYPE double
#include "source_nrm2_c.h"
#undef BASE_TYPE
}

float  gsl_blas_raw_sasum (size_t N, const float X[], size_t incX)
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

double gsl_blas_raw_dasum (size_t N, const double X[], size_t incX)
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

float  gsl_blas_raw_scasum (size_t N, const gsl_complex_packed_array_float X, size_t incX)
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

double gsl_blas_raw_dzasum (size_t N, const gsl_complex_packed_array X, size_t incX)
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


CBLAS_INDEX gsl_blas_raw_isamax (size_t N, const float X[], size_t incX)
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

CBLAS_INDEX gsl_blas_raw_idamax (size_t N, const double X[], size_t incX)
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

CBLAS_INDEX gsl_blas_raw_icamax (size_t N, const gsl_complex_packed_array_float X, size_t incX)
{
  double max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    const double a = fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    if(a > max) {
      max = a;
      result = i;
    }
    i += incX;
  }
  return result;
}

CBLAS_INDEX gsl_blas_raw_izamax (size_t N, const gsl_complex_packed_array X, size_t incX)
{
  double max = 0.0;
  CBLAS_INDEX n;
  CBLAS_INDEX i;
  CBLAS_INDEX result;
  for(n=0; n<N; n++) {
    const double a = fabs(REAL(X, incX, i)) + fabs(IMAG(X, incX, i));
    if(a > max) {
      max = a;
      result = i;
    }
    i += incX;
  }
  return result;
}


void gsl_blas_raw_sswap (size_t N,
                         float X[], size_t incX,
                         float Y[], size_t incY)
{
#define BASE_TYPE float
#include "source_swap_r.h"
#undef BASE_TYPE
/*
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    float tmp = X[i];
    X[i] = Y[j];
    Y[j] = tmp;
    i += incX;
    j += incY;
  }
  */
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
/*
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
  */
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
/*
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    Y[j] = X[i];
    i += incX;
    j += incY;
  }
  */
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
/*
  size_t n;
  size_t i;
  size_t j;
  for(n=0; n<N; n++) {
    REAL(Y, incY, j) = REAL(X, incX, i);
    IMAG(Y, incY, j) = IMAG(X, incX, i);
    i += incX;
    j += incY;
  }
  */
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
/*
  size_t i;

  if(incX == 1 && incY == 1) {
    int m = N % 4;
    for(i=0; i<m; i++) {
      Y[i] += alpha * X[i];
    }
    for(i=m; i<N-3; i += 4) {
      Y[i]   += alpha*X[i];
      Y[i+1] += alpha*X[i+1];
      Y[i+2] += alpha*X[i+2];
      Y[i+3] += alpha*X[i+3];
    }
  }
  else {
    size_t ix = 0;
    size_t iy = 0;
    
    for(i=0; i<N; i++) {
      Y[iy] += alpha*X[ix];
      ix += incX;
      iy += incY;
    }
  }
  */
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
/*
  if (fabs(REAL0(alpha)) + fabs(IMAG0(alpha)) != 0.0 ) {
    size_t i;
    for(i=0; i<N; i++) {
      const float tmpR = REAL0(alpha) * REAL(X, incX, i) - IMAG0(alpha) * IMAG(X, incX, i);
      const float tmpI = REAL0(alpha) * IMAG(X, incX, i) + IMAG0(alpha) * REAL(X, incX, i);
      REAL(Y, incY, i) += tmpR;
      IMAG(Y, incY, i) += tmpI;
    }
  }
  */
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
/*
  float roe   = ( fabs(*a) > fabs(*b) ? *a : *b );
  float scale = fabs(*a) + fabs(*b);
  float r, z;

  if( scale != 0.0 ) {
    float aos = *a/scale;
    float bos = *b/scale;
    r = scale * sqrt(aos*aos + bos*bos);
    r = sign(1.0,roe)*r;
    *c = *a/r;
    *s = *b/r;
    z = 1.0;
    if( fabs(*a) > fabs(*b) ) z = *s;
    if( fabs(*b) >= fabs(*a) && *c != 0.0 ) z = 1.0/(*c);
  }
  else {
   *c = 1.0;
   *s = 0.0;
   r = 0.0;
   z = 0.0;
  }
  
  *a = r;
  *b = z;
  */
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
/*
        SUBROUTINE SROTMG (SD1,SD2,SX1,SY1,SPARAM)
C
C     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
C     THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*
C     SY2)**T.
C     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
C
C     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
C
C       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
C     H=(          )    (          )    (          )    (          )
C       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
C     LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
C     RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
C     VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
C
C     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
C     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
C     OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
C

  const float gam = 4096.0;
  const float gamsq = 1.67772e+07;
  const float rgamsq = 5.96046e-08;

      IF(.NOT. SD1 .LT. ZERO) GO TO 10
C       GO ZERO-H-D-AND-SX1..
          GO TO 60
   10 CONTINUE
C     CASE-SD1-NONNEGATIVE
      SP2=SD2*SY1
      IF(.NOT. SP2 .EQ. ZERO) GO TO 20
          SFLAG=-TWO
          GO TO 260
C     REGULAR-CASE..
   20 CONTINUE
      SP1=SD1*SX1
      SQ2=SP2*SY1
      SQ1=SP1*SX1
C
      IF(.NOT. ABS(SQ1) .GT. ABS(SQ2)) GO TO 40
          SH21=-SY1/SX1
          SH12=SP2/SP1
C
          SU=ONE-SH12*SH21
C
          IF(.NOT. SU .LE. ZERO) GO TO 30
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   30     CONTINUE
               SFLAG=ZERO
               SD1=SD1/SU
               SD2=SD2/SU
               SX1=SX1*SU
C         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF(.NOT. SQ2 .LT. ZERO) GO TO 50
C         GO ZERO-H-D-AND-SX1..
               GO TO 60
   50     CONTINUE
               SFLAG=ONE
               SH11=SP1/SP2
               SH22=SX1/SY1
               SU=ONE+SH11*SH22
               STEMP=SD2/SU
               SD2=SD1/SU
               SD1=STEMP
               SX1=SY1*SU
C         GO SCALE-CHECK
               GO TO 100
C     PROCEDURE..ZERO-H-D-AND-SX1..
   60 CONTINUE
          SFLAG=-ONE
          SH11=ZERO
          SH12=ZERO
          SH21=ZERO
          SH22=ZERO
C
          SD1=ZERO
          SD2=ZERO
          SX1=ZERO
C         RETURN..
          GO TO 220
C     PROCEDURE..FIX-H..
   70 CONTINUE
      IF(.NOT. SFLAG .GE. ZERO) GO TO 90
C
          IF(.NOT. SFLAG .EQ. ZERO) GO TO 80
          SH11=ONE
          SH22=ONE
          SFLAG=-ONE
          GO TO 90
   80     CONTINUE
          SH21=-ONE
          SH12=ONE
          SFLAG=-ONE
   90 CONTINUE
      GO TO IGO,(120,150,180,210)
C     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF(.NOT. SD1 .LE. RGAMSQ) GO TO 130
               IF(SD1 .EQ. ZERO) GO TO 160
               ASSIGN 120 TO IGO
C              FIX-H..
               GO TO 70
  120          CONTINUE
               SD1=SD1*GAM**2
               SX1=SX1/GAM
               SH11=SH11/GAM
               SH12=SH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF(.NOT. SD1 .GE. GAMSQ) GO TO 160
               ASSIGN 150 TO IGO
C              FIX-H..
               GO TO 70
  150          CONTINUE
               SD1=SD1/GAM**2
               SX1=SX1*GAM
               SH11=SH11*GAM
               SH12=SH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF(.NOT. ABS(SD2) .LE. RGAMSQ) GO TO 190
               IF(SD2 .EQ. ZERO) GO TO 220
               ASSIGN 180 TO IGO
C              FIX-H..
               GO TO 70
  180          CONTINUE
               SD2=SD2*GAM**2
               SH21=SH21/GAM
               SH22=SH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF(.NOT. ABS(SD2) .GE. GAMSQ) GO TO 220
               ASSIGN 210 TO IGO
C              FIX-H..
               GO TO 70
  210          CONTINUE
               SD2=SD2/GAM**2
               SH21=SH21*GAM
               SH22=SH22*GAM
          GO TO 200
  220 CONTINUE
          IF(SFLAG)250,230,240
  230     CONTINUE
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               GO TO 260
  240     CONTINUE
               SPARAM(2)=SH11
               SPARAM(5)=SH22
               GO TO 260
  250     CONTINUE
               SPARAM(2)=SH11
               SPARAM(3)=SH21
               SPARAM(4)=SH12
               SPARAM(5)=SH22
  260 CONTINUE
          SPARAM(1)=SFLAG
          RETURN
      END
*/
}

void gsl_blas_raw_drotmg (double d1[], double d2[],
                          double b1[], double b2, double P[])
{
}


void gsl_blas_raw_srot (size_t N,
                        float X[], size_t incX,
                        float Y[], size_t incY,
                        const float c, const float s)
{
#define BASE_TYPE float
#include "source_rot.h"
#undef BASE_TYPE
/*
  size_t i;
  size_t ix = 0;
  size_t iy = 0;
  for(i=0; i<N; i++) {
    float x = X[ix];
    float y = Y[iy];
    X[ix] =  c*x + s*y;
    Y[iy] = -s*x + c*y;
    ix += incX;
    iy += incY;
  }
  */
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
}

void gsl_blas_raw_drotm (size_t N,
                         double X[], size_t incX,
                         double Y[], size_t incY,
                         const double P[])
{
}


void gsl_blas_raw_sscal  (size_t N, float  alpha, float  X[], size_t incX)
{
#include "source_scal_r.h"
/*
  size_t n;
  size_t ix = 0;
  for(n=0; n<N; n++) {
    X[ix] *= alpha;
    ix += incX;
  }
  */
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
/*
  size_t n;
  size_t ix = 0;
  for(n=0; n<N; n++) {
    float tmpR = REAL(X,incX,ix);
    float tmpI = IMAG(X,incX,ix);
    REAL(X,incX,ix) = tmpR*REAL0(alpha) - tmpI*IMAG0(alpha);
    IMAG(X,incX,ix) = tmpR*IMAG0(alpha) + tmpI*REAL0(alpha);
    ix += incX;
  }
  */
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
/*
  size_t n;
  size_t ix = 0;
  for(n=0; n<N; n++) {
    REAL(X,incX,ix) *= alpha;
    IMAG(X,incX,ix) *= alpha;
    ix += incX;
  }
  */
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
/*
  size_t i, j;
  size_t ix, iy;
  size_t lenX, lenY;
  
  if(TransA == CblasNoTrans) {
    lenX = N;
    lenY = M;
  }
  else {
    lenX = M;
    lenY = N;
  }

  iy = 0;
  for(i=0; i<lenY; i++) {
    Y[iy] *= beta;
    iy += incY;
  }

  if(TransA == CblasNoTrans) {
    iy = 0;
    for(i=0; i<lenY; i++) {
      float dot = 0.0;
      ix = 0;
      for(j=0; j<lenX; j++) {
        dot += X[ix]*A[lda*i + j];
        ix += incX;
      }
      Y[iy] += alpha * dot;
      iy += incY;
    }
  }
  else {
    ix = 0;
    for(j=0; j<lenX; j++) {
      float tmp = alpha * X[ix];
      iy = 0;
      for(i=0; i<lenY; i++) {
        Y[iy] += A[lda*j + i] * tmp;
	iy += incY;
      }
      ix += incX;
    }
  }
  */
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
}

void gsl_blas_raw_zgemv (CBLAS_TRANSPOSE TransA,
                         size_t M, size_t N,
                         const gsl_complex_packed alpha,
                         const gsl_complex_packed_array A, int lda,
                         const gsl_complex_packed_array X, size_t incX,
                         const gsl_complex_packed beta,
                         gsl_complex_packed_array Y, size_t incY)
{
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
}

void gsl_blas_raw_dgbmv (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N, size_t KL, size_t KU,
                               double alpha,
                               const double A[], int lda,
                               const double X[], size_t incX,
                               double beta,
                               double Y[], size_t incY)
{
}

void gsl_blas_raw_cgbmv (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N, size_t KL, size_t KU,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float A, int lda,
                               const gsl_complex_packed_array_float X, size_t incX,
                               const gsl_complex_packed_float beta,
                               gsl_complex_packed_array_float Y, size_t incY)
{
}

void gsl_blas_raw_zgbmv (CBLAS_TRANSPOSE TransA,
                               size_t M, size_t N, size_t KL, size_t KU,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, size_t incY)
{
}


/* TRMV */

void gsl_blas_raw_strmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA,
                               CBLAS_DIAG Diag,
                               size_t N,
                               const float A[], int lda,
                               float X[], size_t incX)
{

}

void gsl_blas_raw_dtrmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double A[], int lda,
                               double X[], size_t incX)
{
}

void gsl_blas_raw_ctrmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, size_t incX)
{
}

void gsl_blas_raw_ztrmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, size_t incX)
{
}


/* TBMV */

void gsl_blas_raw_stbmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const float A[], int lda,
                               float X[], size_t incX)
{
}

void gsl_blas_raw_dtbmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const double A[], int lda,
                               double X[], size_t incX)
{
}

void gsl_blas_raw_ctbmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, size_t incX)
{
}

void gsl_blas_raw_ztbmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, size_t incX)
{
}


/* TPMV */

void gsl_blas_raw_stpmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const float Ap[],
                               float X[], size_t incX)
{
}

void gsl_blas_raw_dtpmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double Ap[],
                               double X[], size_t incX)
{
}

void gsl_blas_raw_ctpmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float Ap,
                               gsl_complex_packed_array_float X, size_t incX)
{
}

void gsl_blas_raw_ztpmv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array Ap,
                               gsl_complex_packed_array X, size_t incX)
{
}

/* TRSV */

void gsl_blas_raw_strsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const float A[], int lda,
                               float X[], size_t incX)
{
}

void gsl_blas_raw_dtrsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double A[], int lda,
                               double X[], size_t incX)
{
}

void gsl_blas_raw_ctrsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, size_t incX)
{
}

void gsl_blas_raw_ztrsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, size_t incX)
{
}


/* TBSV */

void gsl_blas_raw_stbsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const float A[], int lda,
                               float X[], size_t incX)
{
}

void gsl_blas_raw_dtbsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const double A[], int lda,
                               double X[], size_t incX)
{
}

void gsl_blas_raw_ctbsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array_float A, int lda,
                               gsl_complex_packed_array_float X, size_t incX)
{
}

void gsl_blas_raw_ztbsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N, size_t K,
                               const gsl_complex_packed_array A, int lda,
                               gsl_complex_packed_array X, size_t incX)
{
}


/* TPSV */

void gsl_blas_raw_stpsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const float Ap[],
                               float X[], size_t incX)
{
}

void gsl_blas_raw_dtpsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const double Ap[],
                               double X[], size_t incX)
{
}

void gsl_blas_raw_ctpsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array_float Ap,
                               gsl_complex_packed_array_float X, size_t incX)
{
}

void gsl_blas_raw_ztpsv (CBLAS_UPLO Uplo,
                               CBLAS_TRANSPOSE TransA, CBLAS_DIAG Diag,
                               size_t N,
                               const gsl_complex_packed_array Ap,
                               gsl_complex_packed_array X, size_t incX)
{
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
/*
  size_t i, j;
  size_t ix, iy, jx, jy;

  iy = 0;
  for(i=0; i<N; i++) {
    Y[iy] *= beta;
    iy += incY;
  }

  if(Uplo == CblasUpper) {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      float tmp1 = alpha * X[jx];
      float tmp2 = 0.0;
      Y[jy] += tmp1*A[lda*j + j];
      ix = jx;
      iy = jy;
      for(i=j+1; i<N; i++) {
        ix += incX;
        iy += incY;
        Y[iy] += tmp1 * A[lda*j + i];
	tmp2  += A[lda*j + i] * X[ix];
      }
      Y[jy] += alpha * tmp2;
      jx += incX;
      jy += incY;
    }
  }
  else {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      float tmp1 = alpha * X[jx];
      float tmp2 = 0.0;
      ix = 0;
      iy = 0;
      for(i=0; i<j; i++) {
        Y[iy] += tmp1 * A[lda*j + i];
	tmp2  += A[lda*j + i] * X[ix];
	ix += incX;
        iy += incY;
      }
      Y[jy] += tmp1*A[lda*j + j] + alpha * tmp2;
      jx += incX;
      jy += incY;
    }
  }
*/
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
  size_t i, j, ell;
  size_t ix, iy, jx, jy;
  size_t kp1 = K + 1;

  iy=0;
  for(i=0; i<N; i++) {
    Y[iy] *= beta;
    iy += incY;
  }
  
  jx = 0;
  jy = 0;
  for(j=0; j<N; j++) {
    float tmp1 = alpha * X[jx];
    float tmp2 = 0.0;
    ix = 0;
    iy = 0;
    ell = kp1 - j;
    for(i=locMAX(0, j-K); i<j; i++) {
      float ajli = A[lda*j + ell+i];
      Y[iy] += tmp1 * ajli;
      tmp2  += ajli * X[ix];
      ix += incX;
      iy += incY;
    }
    Y[jy] += tmp1*A[lda*j + kp1] + alpha*tmp2;
    jx += incX;
    jy += incY;
  }
  
  /* FIXME */
}

void gsl_blas_raw_dsbmv (CBLAS_UPLO Uplo,
                         size_t N, size_t K,
                         double alpha,
                         const double A[], int lda,
                         const double X[], size_t incX,
                         double beta,
                         double Y[], size_t incY)
{
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
/*
  size_t i, j, k;
  size_t ix, iy, jx, jy;
  size_t kk = 0;

  iy=0;
  for(i=0; i<N; i++) {
    Y[iy] *= beta;
    iy += incY;
  }
  
  if(Uplo == CblasUpper) {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      float tmp1 = alpha*X[jx];
      float tmp2 = 0.0;
      ix = 0;
      iy = 0;
      for(k=kk; k<kk+j-1; k++) {
        float apk = Ap[k];
        Y[iy] += tmp1 * apk;
        tmp2  += apk * X[ix];
        ix += incX;
        iy += incY;
      }
      Y[jy] += tmp1*Ap[kk + j - 1] + alpha*tmp2;
      jx += incX;
      jy += incY;
      kk += j;
    }
  }
  else {
    jx = 0;
    jy = 0;
    for(j=0; j<N; j++) {
      float tmp1 = alpha*X[jx];
      float tmp2 = 0.0;
      Y[jy] += tmp1*Ap[kk];
      ix = jx;
      iy = jy;
      for(k=kk+1; k<=kk+N-j; k++) {
        float apk = Ap[k];
        ix += incX;
	iy += incY;
	Y[iy] += tmp1 * apk;
	tmp2  += apk * X[ix];
      }
      Y[jy] += alpha*tmp2;
      jx += incX;
      jy += incY;
      kk += N - j + 1;
    }
  }
*/
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
/*
  size_t i, j;
  size_t ix, jy;
  
  jy = 0;
  for(j=0; j<M; j++) {
    float tmp = alpha * Y[jy];
    ix = 0;
    for(i=0; i<N; i++) {
      A[lda*j + i] += X[ix]*tmp;
      ix += incX;
    }
    jy += incY;
  }
  */
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
/*
  size_t i, j;
  size_t ix, jx;

  if(Uplo == CblasUpper) {
    jx = 0;
    for(j=0; j<N; j++) {
      float tmp = alpha * X[jx];
      ix = jx;
      for(i=j; i<N; i++) {
        A[lda*j + i] += X[ix] * tmp;
        ix += incX;
      }
      jx += incX;
    }
  }
  else {
    jx = 0;
    for(j=0; j<N; j++) {
      float tmp = alpha * X[jx];
      ix = 0;
      for(i=0; i<j; i++) {
        A[lda*j + i] += X[ix]*tmp;
        ix += incX;
      }
      jx += incX;
    }
  }
  */
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
/*
  size_t i, j, k;
  size_t ix, iy, jx, jy, kk;
  
  kk = 0;

  if(Uplo == CblasUpper) {
    jx = 0;
    for(j=0; j<N; j++) {
      float tmp = alpha * X[jx];
      ix = 0;
      for(k=kk; k<kk+j; k++) {
        Ap[k] += X[ix] * tmp;
        ix += incX;
      }
      jx += incX;
      kk += j;
    }
  }
  else {
    jx = 0;
    for(j=0; j<N; j++) {
      float tmp = alpha * X[jx];
      ix = jx;
      for(k=kk; k<=kk+N-j; k++) {
        Ap[k] += X[ix]*tmp;
      }
      jx += incX;
      kk += N - j + 1;
    }
  }
  */
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
/*
  size_t i, j;
  size_t ix, iy, jx, jy;
  
  jx = 0;
  jy = 0;
  
  if(Uplo == CblasUpper) {
    for(j=0; j<N; j++) {
      float tmp1 = alpha * Y[jy];
      float tmp2 = alpha * X[jx];
      ix = jx;
      iy = jy;
      for(i=j; i<N; i++) {
        A[lda*j + i] += X[ix]*tmp1 + Y[iy]*tmp2;
	ix += incX;
	iy += incY;
      }
      jx += incX;
      jy += incY;
    }
  }
  else {
    for(j=0; j<N; j++) {
      float tmp1 = alpha * Y[jy];
      float tmp2 = alpha * X[jx];
      ix = 0;
      iy = 0;
      for(i=0; i<j; i++) {
        A[lda*j + i] += X[ix]*tmp1 + Y[iy]*tmp2;
        ix += incX;
        iy += incY;
      }
      jx += incX;
      jy += incY;
    }
  }
  */
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
}

void gsl_blas_raw_zhemv (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, size_t incY)
{
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
}

void gsl_blas_raw_zhbmv (CBLAS_UPLO Uplo,
                               size_t N, size_t K,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array A, int lda,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, size_t incY)
{
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
}

void gsl_blas_raw_zhpmv (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array Ap,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed beta,
                               gsl_complex_packed_array Y, size_t incY)
{
}


/* GERU */

void gsl_blas_raw_cgeru (size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, size_t incX,
                               const gsl_complex_packed_array_float Y, size_t incY,
                               gsl_complex_packed_array_float A, int lda)
{
}

void gsl_blas_raw_zgeru (size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed_array Y, size_t incY,
                               gsl_complex_packed_array A, int lda)
{
}


/* GERC */

void gsl_blas_raw_cgerc (size_t M, size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, size_t incX,
                               const gsl_complex_packed_array_float Y, size_t incY,
                               gsl_complex_packed_array_float A, int lda)
{
}

void gsl_blas_raw_zgerc (size_t M, size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed_array Y, size_t incY,
                               gsl_complex_packed_array A, int lda)
{
}

/* HER */

void gsl_blas_raw_cher (CBLAS_UPLO Uplo,
                              size_t N,
                              float alpha,
                              const gsl_complex_packed_array_float X, size_t incX,
                              gsl_complex_packed_array_float A, int lda)
{
}

void gsl_blas_raw_zher (CBLAS_UPLO Uplo,
                              size_t N,
                              double alpha,
                              const gsl_complex_packed_array X, size_t incX,
                              gsl_complex_packed_array A, int lda)
{
}


/* HPR */

void gsl_blas_raw_chpr (CBLAS_UPLO Uplo,
                              size_t N,
                              float alpha,
                              const gsl_complex_packed_array_float X, size_t incX,
                              gsl_complex_packed_array_float A)
{
}

void gsl_blas_raw_zhpr (CBLAS_UPLO Uplo,
                              size_t N,
                              double alpha,
                              const gsl_complex_packed_array X, size_t incX,
                              gsl_complex_packed_array A)
{
}


/* HER2 */

void gsl_blas_raw_cher2 (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, size_t incX,
                               const gsl_complex_packed_array_float Y, size_t incY,
                               gsl_complex_packed_array_float A, int lda)
{
}

void gsl_blas_raw_zher2 (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed_array Y, size_t incY,
                               gsl_complex_packed_array A, int lda)
{
}


/* HPR2 */

void gsl_blas_raw_chpr2 (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed_float alpha,
                               const gsl_complex_packed_array_float X, size_t incX,
                               const gsl_complex_packed_array_float Y, size_t incY,
                               gsl_complex_packed_array_float Ap)
{
}

void gsl_blas_raw_zhpr2 (CBLAS_UPLO Uplo,
                               size_t N,
                               const gsl_complex_packed alpha,
                               const gsl_complex_packed_array X, size_t incX,
                               const gsl_complex_packed_array Y, size_t incY,
                               gsl_complex_packed_array Ap)
{
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
