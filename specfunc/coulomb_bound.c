/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_laguerre.h"
#include "gsl_sf_coulomb.h"


/* normalization for hydrogenic wave functions */
static
double
R_norm(const int n, const int l, const double Z)
{
  double A = 2.0*Z/n;
  double term1 = A*A*A /(2.0*n);
  gsl_sf_result ln_a, ln_b;
  gsl_sf_lnfact_impl(n+l, &ln_a);
  gsl_sf_lnfact_impl(n-l-1, &ln_b);
  return sqrt(term1) * exp(0.5*(ln_b.val-ln_a.val));
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
    double norm = R_norm(n, l, Z);
    double rho = A*r;
    double ea = exp(-0.5*rho);
    double pp = gsl_sf_pow_int(rho, l);
    gsl_sf_result lag;
    int stat_lag = gsl_sf_laguerre_n_impl(n-l-1, 2*l+1, rho, &lag);
    int stat_uf;
    double W = norm * ea * pp;
    result->val  = W * lag.val;
    result->err  = fabs(W)*lag.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(0.5*rho+1.0) * fabs(result->val);
    stat_uf = ( result->val == 0.0 ? GSL_EUNDRFLW : GSL_SUCCESS );
    return GSL_ERROR_SELECT_2(stat_lag, stat_uf);
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
