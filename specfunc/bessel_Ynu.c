/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_bessel.h"


int
gsl_sf_bessel_Ynu_impl(double nu, double x, double * result)
{
  if(x <= 0.0 || nu < 0.0) {
    return GSL_EDOM;
  }
  else {
    int N = (int)(nu + 0.5);
    double mu = nu - N;      /* -1/2 <= mu <= 1/2 */
    double Y_mu, Y_mup1, Yp_mu;
    double Y_nu, Yp_nu;

    if(x < 2.0) {
      gsl_sf_bessel_Y_temme(mu, x, &Y_mu, &Y_mup1, &Yp_mu);
    }
    else {
      double J_mu, J_mup1, Jp_mu;
      double P, Q;
      if(mu >= 0.0) {
        gsl_sf_bessel_Jnu_impl(mu,     x, &J_mu);
        gsl_sf_bessel_Jnu_impl(mu+1.0, x, &J_mup1);
      }
      else {
        double J_mup2;
	gsl_sf_bessel_Jnu_impl(mu+1.0, x, &J_mup1);
        gsl_sf_bessel_Jnu_impl(mu+2.0, x, &J_mup2);
	J_mu = 2.0*(mu+1.0)/x * J_mup1 - J_mup2;
      }
      gsl_sf_bessel_JY_steed_CF2(mu, x, &P, &Q);
      Jp_mu = mu/x * J_mu - J_mup1;
      Y_mu  = (P*J_mu - Jp_mu)/Q;
      Yp_mu = (Y_mu*Jp_mu + 2.0/(M_PI*x))/J_mu;
    }
    gsl_sf_bessel_Y_recur(mu, x, N, Y_mu, Yp_mu, &Y_nu, &Yp_nu, (double *)0, (double *)0);
    *result = Y_nu;
    return GSL_SUCCESS;
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
