/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_bessel.h"


int
gsl_sf_bessel_Ynu_impl(double nu, double x, double * result)
{
  int n = rint(nu);

  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(fabs(nu - n) < GSL_MACH_EPS) {
    return gsl_sf_bessel_Yn_impl(n, x, result);
  }
  else {
    double Jnu;
    double Jmnu;
    double t = nu*M_PI;
    
    int stat_nu  = gsl_sf_bessel_Jnu_impl( nu, x, &Jnu);
    int stat_mnu = gsl_sf_bessel_Jnu_impl(-nu, x, &Jmnu);
    /* FIXME: this will fail since we only support nu > 0.0 */
    
    if(stat_nu == GSL_SUCCESS && stat_mnu == GSL_SUCCESS) {
      *result = (Jnu * cos(t) - Jmnu) / sin(t);
      return GSL_SUCCESS;
    }
    else {
      *result = 0.0;
      return GSL_FAILURE;
    }
  }
}


int
gsl_sf_bessel_Ynu_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Ynu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Ynu_e", status);
  }
  return status;
}


double
gsl_sf_bessel_Ynu(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Ynu_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Ynu", status);
  }
  return y;
}
