/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Inu_scaled_impl(double nu, double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x < 0.0 || nu < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x*x < 10.0*(nu+1.0)*GSL_ROOT5_DBL_EPSILON) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 50, GSL_DBL_EPSILON, &b);
    result->val = ex * b;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 100, GSL_DBL_EPSILON, &b);
    result->val = ex * b;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(0.5/(nu*nu + x*x) < GSL_ROOT3_DBL_EPSILON) {
    return gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu, x, result);
  }
  else {
    int N = (int)(nu + 0.5);
    double mu = nu - N;      /* -1/2 <= mu <= 1/2 */ 
    double K_mu, K_mup1, Kp_mu;
    double K_nu, K_nup1, K_num1;
    double I_nu_ratio;
    int stat_Irat;
    int stat_Kmu;
    int n;

    /* obtain K_mu, K_mup1 */
    if(x < 2.0) {
      stat_Kmu = gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }
    else {
      stat_Kmu = gsl_sf_bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }

    /* recurse forward to obtain K_num1, K_nu */
    K_nu   = K_mu;
    K_nup1 = K_mup1;

    for(n=0; n<N; n++) {
      K_num1 = K_nu;
      K_nu   = K_nup1;
      K_nup1 = 2.0*(mu+n+1)/x * K_nu + K_num1;
    }

    /* calculate I_{nu+1}/I_nu */
    stat_Irat = gsl_sf_bessel_I_CF1_ser(nu, x, &I_nu_ratio);

    /* solve for I_nu */
    result->val = 1.0/(x * (K_nup1 + I_nu_ratio * K_nu));
    result->err = GSL_DBL_EPSILON * (0.5*N + 2.0) * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_Kmu, stat_Irat);
  }
}


int
gsl_sf_bessel_Inu_impl(double nu, double x, gsl_sf_result * result)
{
  gsl_sf_result b;
  int stat_I = gsl_sf_bessel_Inu_scaled_impl(nu, x, &b);
  int stat_e = gsl_sf_exp_mult_err_impl(x, fabs(x*GSL_DBL_EPSILON),
                                        b.val, b.err,
					result);
  return GSL_ERROR_SELECT_2(stat_e, stat_I);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Inu_scaled_e(double nu, double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Inu_scaled_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Inu_scaled_e", status);
  }
  return status;
}


int
gsl_sf_bessel_Inu_e(double nu, double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Inu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Inu_e", status);
  }
  return status;
}
