/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_trig.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_clausen.h"


static double aclaus_data[15] = {
  2.142694363766688447e+00,
  0.723324281221257925e-01,
  0.101642475021151164e-02,
  0.3245250328531645e-04,
  0.133315187571472e-05,
  0.6213240591653e-07,
  0.313004135337e-08,
  0.16635723056e-09,
  0.919659293e-11,
  0.52400462e-12,
  0.3058040e-13,
  0.18197e-14,
  0.1100e-15,
  0.68e-17,
  0.4e-18
};
static struct gsl_sf_cheb_series aclaus_cs = {
  aclaus_data,
  14,
  -1, 1,
  (double *)0,
  (double *)0
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_clausen_impl(double x, double *result)
{ 
  const double pi_squared = M_PI*M_PI;
  const double x_cut = M_PI * GSL_SQRT_MACH_EPS;

  double sgn = 1.0;
  int status_red;
  int status = GSL_SUCCESS;

  if(x < 0.0) {
    x   = -x;
    sgn = -1.0;
  }

  /* Argument reduction to [0, 2pi) */
  status_red = gsl_sf_angle_restrict_pos_impl(&x);
  status = status_red;

  /* Further reduction to [0,pi) */
  if(x > M_PI) {
    /* simulated extra precision: 2PI = p0 + p1 */
    const double p0 = 6.28125;
    const double p1 = 0.19353071795864769253e-02;
    x = (p0 - x) + p1;
    sgn = -sgn;
  }

  if(x == 0.0) {
    *result = 0.0;
    return status;
  }
  else if(x < x_cut) {
    *result = x * (1.0 - log(x));
  }
  else {
    const double t = 2.0*(x*x / pi_squared - 0.5);
    const double c = gsl_sf_cheb_eval(&aclaus_cs, t);
    *result = x * (c - log(x));
  }

  *result *= sgn;
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_clausen_e(const double x, double * result)
{
  int status = gsl_sf_clausen_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_clausen_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/

double gsl_sf_clausen(const double x)
{
  double y;
  int status = gsl_sf_clausen_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_clausen", status);
  }
  return y;
}
