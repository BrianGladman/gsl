/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_laguerre.h"
#include "gsl_sf_coulomb.h"


/* normalization for hydrogenic wave functions */
static
int
R_norm(const int n, const int l, const double Z, gsl_sf_result * result)
{
  double A   = 2.0*Z/n;
  double pre = sqrt(A*A*A /(2.0*n));
  gsl_sf_result ln_a, ln_b;
  gsl_sf_result ex;
  int stat_a = gsl_sf_lnfact_impl(n+l, &ln_a);
  int stat_b = gsl_sf_lnfact_impl(n-l-1, &ln_b);
  double diff_val = 0.5*(ln_b.val - ln_a.val);
  double diff_err = 0.5*(ln_b.err + ln_a.err) + GSL_DBL_EPSILON * fabs(diff_val);
  int stat_e = gsl_sf_exp_err_impl(diff_val, diff_err, &ex);
  result->val  = pre * ex.val;
  result->err  = pre * ex.err;
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_ERROR_SELECT_3(stat_e, stat_a, stat_b);
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_hydrogenicR_1_impl(const double Z, const double r, gsl_sf_result * result)
{
  if(Z > 0.0 && r >= 0.0) {
    double A = 2.0*Z;
    double norm = A*sqrt(Z);
    double ea = exp(-Z*r);
    result->val = norm*ea;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) * fabs(Z*r);
    return ( result->val == 0.0 ? GSL_EUNDRFLW : GSL_SUCCESS );
  }
  else {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
}


int
gsl_sf_hydrogenicR_impl(const int n, const int l,
                        const double Z, const double r,
                        gsl_sf_result * result)
{
  if(n < 1 || l > n-1 || Z <= 0.0 || r < 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    double A = 2.0*Z/n;
    gsl_sf_result norm;
    int stat_norm = R_norm(n, l, Z, &norm);
    double rho = A*r;
    double ea = exp(-0.5*rho);
    double pp = gsl_sf_pow_int(rho, l);
    gsl_sf_result lag;
    int stat_lag = gsl_sf_laguerre_n_impl(n-l-1, 2*l+1, rho, &lag);
    int stat_uf;
    double W_val = norm.val * ea * pp;
    double W_err = norm.err * ea * pp;
    W_err += norm.val * ((0.5*rho + 1.0) * GSL_DBL_EPSILON) * ea * pp;
    W_err += norm.val * ea * ((l+1.0) * GSL_DBL_EPSILON) * pp;
    result->val  = W_val * lag.val;
    result->err  = W_val * lag.err + W_err * fabs(lag.val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    stat_uf = ( result->val == 0.0 ? GSL_EUNDRFLW : GSL_SUCCESS );
    return GSL_ERROR_SELECT_3(stat_lag, stat_uf, stat_norm);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_hydrogenicR_1_e(const double Z, const double r, gsl_sf_result * result)
{
  int status = gsl_sf_hydrogenicR_1_impl(Z, r, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hydrogenicR_1_e", status);
  }
  return status;
}


int
gsl_sf_hydrogenicR_e(const int n, const int l, const double Z, const double r, gsl_sf_result * result)
{
  int status = gsl_sf_hydrogenicR_impl(n, l, Z, r, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hydrogenicR_e", status);
  }
  return status;
}
