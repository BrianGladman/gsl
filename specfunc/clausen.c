/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_trig.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_clausen.h"


/* based on MISCFUN CLAUSN(x) by Dr. Allan J. MacLeod [TOMS 757] */

/* CLAUSN(x) = integral 0 to x of (-ln(2*sin(t/2))) dt */

static double aclaus_data[16] = {
  2.1426943637666884470e+00,
  0.7233242812212579245e-01,
  0.101642475021151164e-02,
  0.3245250328531645e-04,
  0.133315187571472e-05,
  0.6213240591653e-07,
  0.313004135337e-08,
  0.16635723056e-09,
  0.919659293e-11,
  0.52400462e-12,
  0.3058040e-13,
  0.181969e-14,
  0.11004e-15,
  0.675e-17,
  0.42e-18,
  0.3e-19
};
static struct gsl_sf_ChebSeries aclaus_cs = {
  aclaus_data,
  15,
  -1, 1,
  (double *)0,
  (double *)0
};


int gsl_sf_clausen_impl(double x, double *result)
{
  /* simulated extra precision constants: 2PI = 2PI_a + 2PI_b */
  const double two_pi_a = 6.28125;
  const double two_pi_b = 0.19353071795864769253e-02;

  const double PISQ     = 9.8696044010893586188;
  double xlo = M_PI * GSL_SQRT_MACH_EPS;
  double sgn = 1.0;
  int status_red;

  if(x < 0.0) {
    x   = -x;
    sgn = -1.0;
  }

  /* Argument reduction to [0, 2pi) */
  status_red = gsl_sf_angle_restrict_pos(&x, 100.0*GSL_MACH_EPS);
  if(status_red != GSL_SUCCESS) {
    *result = 0.0;
    return status_red;
  }

  /* Further reduction to [0,pi) */
  if(x > M_PI) {
    x = (two_pi_a - x) + two_pi_b;
    sgn = -sgn;
  }

  /* Set result to zero if X multiple of PI */
  if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }

  if(x < xlo) {
    *result = x * (1.0 - log(x));
  }
  else { /* xlo <= x < PI */
    double t = 2.0*(x*x / PISQ - 0.5);
    if(t > 1.0) t = 1.0;
    *result = x * gsl_sf_cheb_eval(t, &aclaus_cs) - x * log(x);
  }
  
  *result *= sgn;
  return GSL_SUCCESS;
}


int gsl_sf_clausen_e(double x, double * result)
{
  int status = gsl_sf_clausen_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_clausen_e", status);
  }
  return status;
}

double gsl_sf_clausen(double x)
{
  double y;
  int status = gsl_sf_clausen_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_clausen_e", status);
  }
  return y;
}
