/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"

#define locMin(a,b)  ((a) < (b) ? (a) : (b))
#define locMax(a,b)  ((a) > (b) ? (a) : (b))


int
gsl_sf_bessel_Knu_scaled_impl(const double nu, const double x, double * result)
{
  if(x <= 0.0 || nu < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    int N = (int)(nu + 0.5);
    double mu = nu - N;      /* -1/2 <= mu <= 1/2 */
    double K_mu, K_mup1, Kp_mu;
    double K_nu, Kp_nu;

    if(x < 2.0) {
      gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }
    else {
      gsl_sf_bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }
    gsl_sf_bessel_K_recur(mu, x, N, K_mu, Kp_mu, &K_nu, &Kp_nu, (double *)0, (double *)0);
    *result = K_nu;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_bessel_Knu_scaled_e(const double nu, const double x, double * result)
{
  int status = gsl_sf_bessel_Knu_scaled_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Knu_scaled_e", status);
  }
  return status;
}


double
gsl_sf_bessel_Knu_scaled(const double nu, const double x)
{
  double y;
  int status = gsl_sf_bessel_Knu_scaled_impl(nu, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Knu_scaled", status);
  }
  return y;
}
