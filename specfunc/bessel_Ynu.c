/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_bessel.h"


/* Perform forward recurrence for Y_nu(x) and Y'_nu(x)
 *
 *        Y_{nu+1} =  nu/x Y_nu - Y'_nu
 *       Y'_{nu+1} = -(nu+1)/x Y_{nu+1} + Y_nu
 */
static
int
bessel_Y_recur(const double nu_min, const double x, const int kmax,
               const double Y_start, const double Yp_start,
	       double * Y_end, double * Yp_end)
{
  double x_inv = 1.0/x;
  double nu = nu_min;
  double Y_nu  = Y_start;
  double Yp_nu = Yp_start;
  int k;

  for(k=1; k<=kmax; k++) {
    double nuox = nu*x_inv;
    double Y_nu_save = Y_nu;
    Y_nu  = -Yp_nu + nuox * Y_nu;
    Yp_nu = Y_nu_save - (nuox+x_inv) * Y_nu;
    nu += 1.0;
  }
  *Y_end  = Y_nu;
  *Yp_end = Yp_nu;
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Ynu_impl(double nu, double x, double * result,
                       const gsl_prec_t goal, const unsigned int err_bits)
{
  if(x <= 0.0 || nu < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    int N = (int)(nu + 0.5);
    double mu = nu - N;      /* -1/2 <= mu <= 1/2 */
    double Y_mu, Y_mup1, Yp_mu;
    double Y_nu, Yp_nu;
    int stat_mu;
    int stat_0, stat_1;

    /* Determine Y_mu, Y'_mu for |mu| < 1/2
     */
    if(x < 2.0) {
      stat_mu = gsl_sf_bessel_Y_temme(mu, x, &Y_mu, &Y_mup1, &Yp_mu);
      stat_0 = GSL_SUCCESS;
      stat_1 = GSL_SUCCESS;
    }
    else {
      double J_mu, J_mup1, Jp_mu;
      double P, Q;
      if(mu >= 0.0) {
        stat_0 = gsl_sf_bessel_Jnu_impl(mu,     x, &J_mu, goal, err_bits);
        stat_1 = gsl_sf_bessel_Jnu_impl(mu+1.0, x, &J_mup1, goal, err_bits);
      }
      else {
        double J_mup2;
	stat_0 = gsl_sf_bessel_Jnu_impl(mu+1.0, x, &J_mup1, goal, err_bits);
        stat_1 = gsl_sf_bessel_Jnu_impl(mu+2.0, x, &J_mup2, goal, err_bits);
	J_mu = 2.0*(mu+1.0)/x * J_mup1 - J_mup2;
      }
      stat_mu = gsl_sf_bessel_JY_steed_CF2(mu, x, &P, &Q);
      Jp_mu = mu/x * J_mu - J_mup1;
      Y_mu  = (P*J_mu - Jp_mu)/Q;
      Yp_mu = (Y_mu*Jp_mu + 2.0/(M_PI*x))/J_mu;
    }

    bessel_Y_recur(mu, x, N, Y_mu, Yp_mu, &Y_nu, &Yp_nu);

    *result = Y_nu;
    return GSL_ERROR_SELECT_3(stat_mu, stat_0, stat_1);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Ynu_e(const double nu, const double x, double * result,
                    const gsl_prec_t goal, const unsigned int err_bits)
{
  int status = gsl_sf_bessel_Ynu_impl(nu, x, result, goal, err_bits);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Ynu_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/

double
gsl_sf_bessel_Ynu(const double nu, const double x,
                  const gsl_prec_t goal, const unsigned int err_bits)
{
  double y;
  int status = gsl_sf_bessel_Ynu_impl(nu, x, &y, goal, err_bits);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Ynu", status);
  }
  return y;
}
