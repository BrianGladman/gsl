/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_bessel.h"

#define locMin(a,b)  ((a) < (b) ? (a) : (b))


int
gsl_sf_bessel_Inu_impl(double nu, double x, double * result)
{
  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_MACH_EPS) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 4, &b);
    *result = ex * b;
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 100, result);
    *result = ex * b;
    return stat;
  }
  else if(locMin( 0.29/(nu*nu), 0.5/(nu*nu + x*x) ) < GSL_ROOT3_MACH_EPS) {
    return gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu, x, result);
  }
  else {
    /* FIXME */
  }
}


int
gsl_sf_bessel_Inu_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Inu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Inu_e", status);
  }
  return status;
}


double
gsl_sf_bessel_Inu(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Inu_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Inu", status);
  }
  return y;
}
