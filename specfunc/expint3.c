/* Author: G. Jungman
 * RCS: $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_expint.h"


static double expint3_data[24] = {
  1.269198414221126014,
 -0.248846446384140982,
  0.80526220717231041e-01,
 -0.25772733251968330e-01,
  0.7599878873073774e-02,
 -0.2030695581940405e-02,
  0.490834586699330e-03,
 -0.107682239142021e-03,
  0.21551726264290e-04,
 -0.3956705137384e-05,
  0.6699240933896e-06,
 -0.105132180807e-06,
  0.15362580199e-07,
 -0.20990960364e-08,
  0.2692109538e-09,
 -0.325195242e-10,
  0.37114816e-11,
 -0.4013652e-12,
  0.412334e-13,
 -0.40338e-14,
  0.3766e-15,
 -0.336e-16,
  0.29e-17,
 -0.2e-18
};
static gsl_sf_cheb_series expint3_cs = {
  expint3_data,
  23,
  -1.0, 1.0,
  (double *)0,
  (double *)0,
  15
};

static double expint3a_data[23] = {
   1.9270464955068273729,
  -0.349293565204813805e-01,
   0.14503383718983009e-02,
  -0.8925336718327903e-04,
   0.70542392191184e-05,
  -0.6671727454761e-06,
   0.724267589982e-07,
  -0.87825825606e-08,
   0.11672234428e-08,
  -0.1676631281e-09,
   0.257550158e-10,
  -0.41957888e-11,
   0.7201041e-12,
  -0.1294906e-12,
   0.24287e-13,
  -0.47331e-14,
   0.95531e-15,
  -0.1991e-15,
   0.428e-16,
  -0.94e-17,
   0.21e-17,
  -0.5e-18,
   0.1e-18
};
static gsl_sf_cheb_series expint3a_cs = {
  expint3a_data,
  22,
  -1.0, 1.0,
  (double *)0,
  (double *)0,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_expint_3_impl(const double x, double * result)
{
  const double val_infinity = 0.892979511569249211;

  if(x < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.6*GSL_ROOT3_MACH_EPS) {
    *result = x;
    return GSL_SUCCESS;
  }
  else if(x <= 2.0) {
    const double t = x*x*x/4.0 - 1.0;
    const double c = gsl_sf_cheb_eval(&expint3_cs, t);
    *result = x * c;
    return GSL_SUCCESS;
  }
  else if(x < pow(-GSL_LOG_MACH_EPS, 1.0/3.0)) {
    const double t = 16.0/(x*x*x) - 1.0;
    const double c = gsl_sf_cheb_eval(&expint3a_cs, t);
    *result = val_infinity - c * exp(-x*x*x)/(3.0*x*x);
    return GSL_SUCCESS;
  }
  else {
    *result = val_infinity;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_expint_3_e(double x, double * result)
{
  int status = gsl_sf_expint_3_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_expint_3_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/

double gsl_sf_expint_3(double x)
{
  double y;
  int status = gsl_sf_expint_3_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_expint_3", status);
  }
  return y;
}
