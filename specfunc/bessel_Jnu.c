/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_bessel.h"


int
gsl_sf_bessel_Jnu_impl(double nu, double x, double * result)
{
  const double nu_cut = 20.0;

  if(x < 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_MACH_EPS) {
    return gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, -1, 4, result);
  }
  else if(x*x < 10.0*(nu+1.0)) {
    return gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, -1, 100, result);
  }
  else if(nu > nu_cut) {
    return gsl_sf_bessel_Jnu_asymp_Olver_impl(nu, x, result);
  }
  else {
    /* Evaluate at large enough nu and apply backward recurrence */

    int n;
    int steps = ceil(nu_cut - nu) + 1;
    double Jnp1;
    double Jn;
    double Jnm1;
    
    gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps + 1.0, x, &Jnp1);
    gsl_sf_bessel_Jnu_asymp_Olver_impl(nu + steps      , x, &Jn);
    
    for(n=steps; n>0; n--) {
      Jnm1 = 2.0*(nu+n)/x * Jn - Jnp1;
      Jnp1 = Jn;
      Jn   = Jnm1;
    }
    
    *result = Jnm1;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_bessel_Jnu_e(double nu, double x, double * result)
{
  int status = gsl_sf_bessel_Jnu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Jnu_e", status);
  }
  return status;
}


double
gsl_sf_bessel_Jnu(double nu, double x)
{
  double y;
  int status = gsl_sf_bessel_Jnu_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Jnu", status);
  }
  return y;
}
