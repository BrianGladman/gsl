/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_expint.h"


static double atanint_data[21] = {
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
};
static gsl_sf_cheb_series atanint_cs = {
  atanint_data,
  20,
  -1, 1,
  (double *)0,
  (double *)0,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_atanint_impl(const double x, gsl_sf_result * result)
{
  const double ax  = fabs(x);
  const double sgn = GSL_SIGN(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(ax == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax < 0.5*GSL_SQRT_DBL_EPSILON) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(ax <= 1.0) {
    const double t = 2.0 * (x*x - 0.5);
    gsl_sf_result result_c;
    gsl_sf_cheb_eval_impl(&atanint_cs, t, &result_c);
    result->val  = x * result_c.val;
    result->err  = x * result_c.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(ax < 1.0/GSL_SQRT_DBL_EPSILON) {
    const double t = 2.0 * (1.0/(x*x) - 0.5);
    gsl_sf_result result_c;
    gsl_sf_cheb_eval_impl(&atanint_cs, t, &result_c);
    result->val  = sgn * (0.5*M_PI*log(ax) + result_c.val/ax);
    result->err  = result_c.err/ax + fabs(result->val*GSL_DBL_EPSILON);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = sgn * 0.5*M_PI*log(ax);
    result->err = 2.0 * fabs(result->val * GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_atanint_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_atanint_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_atanint_e", status);
  }
  return status;
}
