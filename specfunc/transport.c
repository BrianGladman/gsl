/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_transport.h"



/* based on TRANS0n(X) from MISCFUN   by ALLAN J. MACLEOD  [TOMS 757] */

/* transport(n,x) := Integral[ t^n e^t /(e^t - 1)^2, {t,0,x}]  */


static double tran2_data[20] = {
   1.67176044643453850301e+00,
  -0.14773535994679448986e+00,
   0.1482138199469363384e-01,
  -0.141953303263056126e-02,
   0.13065413244157083e-03,
  -0.1171557958675790e-04,
   0.103334984457557e-05,
  -0.9019113042227e-07,
   0.781771698331e-08,
  -0.67445656840e-09,
   0.5799463945e-10,
  -0.497476185e-11,
   0.42596097e-12,
  -0.3642189e-13,
   0.311086e-14,
  -0.26547e-15,
   0.2264e-16,
  -0.193e-17,
   0.16e-18,
  -0.1e-19
};
static struct gsl_sf_cheb_series tran2_cs = {
  tran2_data,
  19,
  -1, 1,
  (double *)0,
  (double *)0
};

static double tran3_data[20] = {
   0.76201254324387200657e+00,
  -0.10567438770505853250e+00,
   0.1197780848196578097e-01,
  -0.121440152036983073e-02,
   0.11550997693928547e-03,
  -0.1058159921244229e-04,
   0.94746633853018e-06,
  -0.8362212128581e-07,
   0.731090992775e-08,
  -0.63505947788e-09,
   0.5491182819e-10,
  -0.473213954e-11,
   0.40676948e-12,
  -0.3489706e-13,
   0.298923e-14,
  -0.25574e-15,
   0.2186e-16,
  -0.187e-17,
   0.16e-18,
  -0.1e-19
};
static struct gsl_sf_cheb_series tran3_cs = {
  tran3_data,
  19,
  -1, 1,
  (double *)0,
  (double *)0
};


static
double
transport_sumexp(int numexp, int order, double t, double x)
{
  double rk = (double)numexp;
  double sumexp = 0.0;
  double xk;
  double xk1;
  int k;
  for(k=1; k<=numexp; k++) {
    int j;
    double sum2 = 1.0;
    xk = 1.0/(rk*x);
    xk1 = 1.0;
    for(j=1; j<=order; j++) {
      sum2 = sum2*xk1*xk + 1.0;
      xk1 += 1.0;
    }
    sumexp = sumexp*t + sum2;
    rk -= 1.0;
  }
  return sumexp;
}


int
gsl_sf_transport_2_impl(const double x, double * result)
{
/*
C  MACHINE-DEPENDENT CONSTANTS:
C
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN02 = x to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large x contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN02 = VALINF  -  x**2 exp(-x)
C                    The recommended value is 2/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.
C
C     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT
*/

  const double valinf = 0.32898681336964528729e+01;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 3.0*GSL_SQRT_MACH_EPS) {
    *result = x;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double t = (x*x/8.0 - 0.5) - 0.5;
    *result = x * gsl_sf_cheb_eval(&tran2_cs, t);
    return GSL_SUCCESS;
  }
  else if(x < -GSL_LOG_MACH_EPS) {
    int    numexp = (int)((-GSL_LOG_MACH_EPS)/x) + 1;
    double sumexp = transport_sumexp(numexp, 2, exp(-x), x);
    double t = 2.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = valinf;
    }
    else {
      *result = valinf - exp(t);
    }
    return GSL_SUCCESS;
  }
  else if(x < 2.0/GSL_MACH_EPS) {
    int    numexp = 1;
    double sumexp = transport_sumexp(numexp, 2, 1.0, x);
    double t = 2.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = valinf;
    }
    else {
      *result = valinf - exp(t);
    }
    return GSL_SUCCESS;
  }
  else {
    double t = 2.0 * log(x) - x;
    if(t < GSL_LOG_MACH_EPS) {
      *result = valinf;
    }
    else {
      *result = valinf - exp(t);
    }
    return GSL_SUCCESS;
  }
}


int
gsl_sf_transport_3_impl(const double x, double * result)
{
/*
C    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
C                       The recommended value is such that
C                             ATRAN(NTERMS) < EPS/100
C
C    XLOW2 - DOUBLE PRECISION - The value below which TRAN03 = 0.0 to machine
C                    precision. The recommended value is
C                          square root of (2*XMIN)
C
C    XLOW1 - DOUBLE PRECISION - The value below which TRAN03 = X**2/2 to
C                   machine precision. The recommended value is
C                             sqrt(8*EPSNEG)
C
C    XHIGH1 - DOUBLE PRECISION - The value above which the exponential series for
C                    large X contains only one term. The recommended value
C                    is        - ln(EPS).
C
C    XHIGH2 - DOUBLE PRECISION - The value above which 
C                       TRAN03 = VALINF  -  X**3 exp(-X)
C                    The recommended value is 3/EPS
C
C    XHIGH3 - DOUBLE PRECISION - The value of ln(EPSNEG). Used to prevent overflow
C                    for large x.


      DATA NUMJN,RNUMJN/ 3 , 3.0 D 0 /

      DATA NTERMS/17/
      DATA XLOW1,XLOW2/2.98023224D-8,2.10953733D-154/
      DATA XHIGH1,XHIGH3/36.04365668D0,-36.73680056D0/
      DATA XHIGH2/1.35107988D16/
*/
  
  const double valinf = 0.72123414189575657124e+01;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0*GSL_SQRT_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(x < 3.0*GSL_SQRT_MACH_EPS) {
    *result = 0.5*x*x;
    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    double x2 = x*x;
    double t = (x2/8.0 - 0.5) - 0.5;
    *result = x2 * gsl_sf_cheb_eval(&tran3_cs, t);
    return GSL_SUCCESS;
  }
  else if(x < -GSL_LOG_MACH_EPS) {
    int    numexp = (int)((-GSL_LOG_MACH_EPS)/x) + 1;
    double sumexp = transport_sumexp(numexp, 3, exp(-x), x);
    double t = 2.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = valinf;
    }
    else {
      *result = valinf - exp(t);
    }
    return GSL_SUCCESS;
  }
  else if(x < 3.0/GSL_MACH_EPS) {
    int    numexp = 1;
    double sumexp = transport_sumexp(numexp, 3, 1.0, x);
    double t = 2.0 * log(x) - x + log(sumexp);
    if(t < GSL_LOG_MACH_EPS) {
      *result = valinf;
    }
    else {
      *result = valinf - exp(t);
    }
    return GSL_SUCCESS;
  }
  else {
    double t = 2.0 * log(x) - x;
    if(t < GSL_LOG_MACH_EPS) {
      *result = valinf;
    }
    else {
      *result = valinf - exp(t);
    }
    return GSL_SUCCESS;
  }
}


int gsl_sf_transport_2_e(const double x, double * result)
{
  int status = gsl_sf_transport_2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_transport_2_e", status);
  }
  return status;
}

int gsl_sf_transport_3_e(const double x, double * result)
{
  int status = gsl_sf_transport_3_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_transport_3_e", status);
  }
  return status;
}


double gsl_sf_transport_2(const double x)
{
  double y;
  int status = gsl_sf_transport_2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_transport_2", status);
  }
  return y;
}

double gsl_sf_transport_3(const double x)
{
  double y;
  int status = gsl_sf_transport_3_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_transport_3", status);
  }
  return y;
}
