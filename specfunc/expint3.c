/* Author: G. Jungman
 * RCS: $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_expint.h"

/* based on MISCFUN EXP3(X) by Allan J. Macleod  [TOMS 757] */

/* expint_3(x) = integral 0 to x  (exp(-t*t*t)) dt */

static double aexp3_data[25] = {
  1.26919841422112601434,
 -0.24884644638414098226,
  0.8052622071723104125e-01,
 -0.2577273325196832934e-01,
  0.759987887307377429e-02,
 -0.203069558194040510e-02,
  0.49083458669932917e-03,
 -0.10768223914202077e-03,
  0.2155172626428984e-04,
 -0.395670513738429e-05,
  0.66992409338956e-06,
 -0.10513218080703e-06,
  0.1536258019825e-07,
 -0.209909603636e-08,
  0.26921095381e-09,
 -0.3251952422e-10,
  0.371148157e-11,
 -0.40136518e-12,
  0.4123346e-13,
 -0.403375e-14,
  0.37658e-15,
 -0.3362e-16,
  0.288e-17,
 -0.24e-18,
  0.2e-19
};
static struct gsl_sf_cheb_series aexp3_cs = {
  aexp3_data,
  24,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};

static double aexp3a_data[25] = {
   1.92704649550682737293,
  -0.3492935652048138054e-01,
   0.145033837189830093e-02,
  -0.8925336718327903e-04,
   0.705423921911838e-05,
  -0.66717274547611e-06,
   0.7242675899824e-07,
  -0.878258256056e-08,
   0.116722344278e-08,
  -0.16766312812e-09,
   0.2575501577e-10,
  -0.419578881e-11,
   0.72010412e-12,
  -0.12949055e-12,
   0.2428703e-13,
  -0.473311e-14,
   0.95531e-15,
  -0.19914e-15,
   0.4277e-16,
  -0.944e-17,
   0.214e-17,
  -0.50e-18,
   0.12e-18,
  -0.3e-19,
   0.1e-19
};
static struct gsl_sf_cheb_series aexp3a_cs = {
  aexp3a_data,
  24,
  -1.0, 1.0,
  (double *)0,
  (double *)0
};
      
int gsl_sf_expint_3_impl(const double x, double * result)
{
  const double FUNINF = 0.89297951156924921122;
  const double xlo = pow(4.0*GSL_MACH_EPS,  1.0/3.0);
  const double xup = pow(-GSL_LOG_MACH_EPS, 1.0/3.0);

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < xlo) {
    *result = x;
    return GSL_SUCCESS;
  }
  else if(x <= 2.0) {
    double t = x*x*x/4.0 - 1.0;
    *result = x * gsl_sf_cheb_eval(t, &aexp3_cs);
    return GSL_SUCCESS;
  }
  else if(x < xup) {
    double t = 16.0/(x*x*x) - 1.0;
    double r = gsl_sf_cheb_eval(t, &aexp3a_cs);
    r *= exp(-x*x*x)/(3.0*x*x);
    *result = FUNINF - r;
    return GSL_SUCCESS;
  }
  else {
    *result = FUNINF;
    return GSL_SUCCESS;
  }
}
