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


/* Perform backward recurrence for I_nu(x) and I'_nu(x)
 *
 *        I_{nu-1} = nu/x I_nu + I'_nu
 *       I'_{nu-1} = (nu-1)/x I_{nu-1} + I_nu
 */
static
int
bessel_I_recur(const double nu_min, const double x, const int kmax,
               const double I_start, const double Ip_start,
	       double * I_end, double * Ip_end
	       )
{
  double x_inv  = 1.0/x;
  double nu_max = nu_min + kmax;
  double I_nu  = I_start;
  double Ip_nu = Ip_start;
  double nu = nu_max;
  int k;

  for(k=kmax-1; k>=0; k--) {
    double nuox = nu*x_inv;
    double I_nu_save = I_nu;
    I_nu  = nuox*I_nu + Ip_nu;
    Ip_nu = (nuox - x_inv)*I_nu + I_nu_save;
    nu -= 1.0;
  }
  *I_end  = I_nu;
  *Ip_end = Ip_nu;
  return GSL_SUCCESS;
}


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
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(x*x < 10.0*(nu+1.0)) {
    double b;
    double ex = exp(-x);
    int stat = gsl_sf_bessel_Inu_Jnu_taylor_impl(nu, x, 1, 100, GSL_DBL_EPSILON, &b);
    result->val = ex * b;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(0.5/(nu*nu + x*x) < GSL_ROOT3_DBL_EPSILON) {
    return gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nu, x, result);
  }
  else {
    int N = floor(nu);
    int M = ceil(sqrt(0.5/GSL_ROOT3_DBL_EPSILON - x*x + 1.0)) + 1;
    double nu_frac = nu - N;
    double mu = M + nu_frac;
    double I_mu, Ip_mu, I_mup1;
    double I_nu, Ip_nu;
    gsl_sf_result r_I_mu;
    gsl_sf_result r_I_mup1;
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl(mu,     x, &r_I_mu);
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl(mu+1.0, x, &r_I_mup1);
    I_mu   = r_I_mu.val;
    I_mup1 = r_I_mup1.val;
    Ip_mu = mu/x * I_mu + I_mup1;
    bessel_I_recur(nu, x, M-N, I_mu, Ip_mu, &I_nu, &Ip_nu);
    result->val = I_nu;
    result->err = I_nu * (GSL_DBL_EPSILON + r_I_mu.err/r_I_mu.val + r_I_mup1.err/r_I_mup1.val);
    return GSL_SUCCESS;
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
