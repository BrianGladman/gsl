/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "bessel.h"
#include "bessel_temme.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_chebyshev.h"
#include "gsl_sf_bessel.h"



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Knu_scaled_impl(const double nu, const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0 || nu < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    int N = (int)(nu + 0.5);
    double mu = nu - N;      /* -1/2 <= mu <= 1/2 */
    double K_mu, K_mup1, Kp_mu;
    double K_nu, K_nup1, K_num1;
    int n;

    if(x < 2.0) {
      gsl_sf_bessel_K_scaled_temme(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }
    else {
      gsl_sf_bessel_K_scaled_steed_temme_CF2(mu, x, &K_mu, &K_mup1, &Kp_mu);
    }

    /* recurse forward to obtain K_num1, K_nu */
    K_nu   = K_mu;
    K_nup1 = K_mup1;

    for(n=0; n<N; n++) {
      K_num1 = K_nu;
      K_nu   = K_nup1;
      K_nup1 = 2.0*(mu+n+1)/x * K_nu + K_num1;
    }

    result->val = K_nu;
    result->err = 2.0 * GSL_DBL_EPSILON * (N + 4.0) * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_bessel_Knu_impl(const double nu, const double x, gsl_sf_result * result)
{
  gsl_sf_result b;
  int stat_K = gsl_sf_bessel_Knu_scaled_impl(nu, x, &b);
  int stat_e = gsl_sf_exp_mult_err_impl(-x, 0.0, b.val, b.err, result);
  return GSL_ERROR_SELECT_2(stat_e, stat_K);
}


int
gsl_sf_bessel_lnKnu_impl(const double nu, const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0 || nu < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(nu == 0.0) {
    gsl_sf_result K_scaled;
    /* This cannot underflow, and
     * it will not throw GSL_EDOM
     * since that is already checked.
     */
    gsl_sf_bessel_K0_scaled_impl(x, &K_scaled);
    result->val  = -x + log(fabs(K_scaled.val));
    result->err  = GSL_DBL_EPSILON * fabs(x) + fabs(K_scaled.err/K_scaled.val);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 2.0 && nu > 1.0) {
    /* Make use of the inequality
     * Knu(x) <= 1/2 (2/x)^nu Gamma(nu),
     * which follows from the integral representation
     * [Abramowitz+Stegun, 9.6.23 (2)]. With this
     * we decide whether or not there is an overflow
     * problem because x is small.
     */
    double ln_bound;
    gsl_sf_result lg_nu;
    gsl_sf_lngamma_impl(nu, &lg_nu);
    ln_bound = -M_LN2 - nu*log(0.5*x) + lg_nu.val;
    if(ln_bound > GSL_LOG_DBL_MAX - 20.0) {
      /* x must be very small or nu very large (or both).
       */
      double xi  = 0.25*x*x;
      double sum = 1.0 - xi/(nu-1.0);
      if(nu > 2.0) sum +=  (xi/(nu-1.0)) * (xi/(nu-2.0));
      result->val  = ln_bound + log(sum);
      result->err  = lg_nu.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    /* can drop-through here */
  }


  {
    /* We passed the above tests, so no problem.
     * Evaluate as usual. Note the possible drop-through
     * in the above code!
     */
    gsl_sf_result K_scaled;
    gsl_sf_bessel_Knu_scaled_impl(nu, x, &K_scaled);
    result->val  = -x + log(fabs(K_scaled.val));
    result->err  = GSL_DBL_EPSILON * fabs(x) + fabs(K_scaled.err/K_scaled.val);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Knu_scaled_e(const double nu, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Knu_scaled_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Knu_scaled_e", status);
  }
  return status;
}
int
gsl_sf_bessel_Knu_e(const double nu, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Knu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Knu_e", status);
  }
  return status;
}

int
gsl_sf_bessel_lnKnu_e(const double nu, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_lnKnu_impl(nu, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_lnKnu_e", status);
  }
  return status;
}
