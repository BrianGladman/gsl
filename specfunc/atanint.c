/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_expint.h"


static double atanint_data[23] = {
  1.91040361296235937512,
 -0.4176351437656746940e-01,
  0.275392550786367434e-02,
 -0.25051809526248881e-03,
  0.2666981285121171e-04,
 -0.311890514107001e-05,
  0.38833853132249e-06,
 -0.5057274584964e-07,
  0.681225282949e-08,
 -0.94212561654e-09,
  0.13307878816e-09,
 -0.1912678075e-10,
  0.278912620e-11,
 -0.41174820e-12,
  0.6142987e-13,
 -0.924929e-14,
  0.140387e-14,
 -0.21460e-15,
  0.3301e-16,
 -0.511e-17,
  0.79e-18,
 -0.12e-18,
  0.2e-19
};
static struct gsl_sf_cheb_series atanint_cs = {
  atanint_data,
  22,
  -1, 1,
  (double *)0,
  (double *)0
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_atanint_impl(double xin, double * result)
{
/*
C   NTERMS - INTEGER - The no. of terms of the array ATNINTT.
C                      The recommended value is such that
C                          ATNINA(NTERMS) < EPS/100   
C
C   XLOW - DOUBLE PRECISION - A bound below which ATNINT(x) = x to machine
C                 precision. The recommended value is
C                     sqrt(EPSNEG/2).
C 
C   XUPPER - DOUBLE PRECISION - A bound on x, above which, to machine precision 
C                   ATNINT(x) = (pi/2)ln x
C                   The recommended value is 1/EPS.

      DATA NTERMS/19/
      DATA XLOW,XUPPER/7.4505806D-9,4.5036D15/
*/
  const double x   = fabs(xin);
  const double sgn = GSL_SIGN(xin);

  if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 0.5*GSL_SQRT_MACH_EPS) {
    *result = xin;
    return GSL_SUCCESS;
  }
  else if(x <= 1.0) {
    double t = 2.0 * (x*x - 0.5);
    *result = xin * gsl_sf_cheb_eval(&atanint_cs, t);
    return GSL_SUCCESS;
  }
  else if(x < 1.0/GSL_SQRT_MACH_EPS) {
    double t = 2.0 * (1.0/(x*x) - 0.5);
    *result = sgn * (0.5*M_PI*log(x) + gsl_sf_cheb_eval(&atanint_cs, t)/x);
    return GSL_SUCCESS;
  }
  else {
    *result = sgn * 0.5*M_PI*log(x);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_atanint_e(double x, double * result)
{
  int status = gsl_sf_atanint_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_atanint_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/


double
gsl_sf_atanint(double x)
{
  double y;
  int status = gsl_sf_atanint_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_atanint", status);
  }
  return y;
}
