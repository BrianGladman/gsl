/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_log.h"
#include "gsl_sf_psi.h"
#include "gsl_sf_gamma.h"

#define locMAX(a,b) ((a) > (b) ? (a) : (b))
#define locMIN(a,b) ((a) < (b) ? (a) : (b))


int
gsl_sf_lnbeta_impl(const double x, const double y, double * result)
{
  if(x <= 0.0 || y <= 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    const double max = locMAX(x,y);
    const double min = locMIN(x,y);
    const double rat = min/max;

    if(rat < 0.1) {
      double lnopr;
      double lnpow;
      double lnpre;
      double gsx, gsy, gsxy;
      gsl_sf_gammastar_impl(x, &gsx);
      gsl_sf_gammastar_impl(y, &gsy);
      gsl_sf_gammastar_impl(x+y, &gsxy);
      gsl_sf_log_1plusx_impl(rat, &lnopr);
      lnpre = log(gsx*gsy/gsxy * M_SQRT2*M_SQRTPI);
      lnpow = min*log(rat) - 0.5*log(min) - (x+y-0.5)*lnopr;
      *result = lnpre + lnpow;
    }
    else {
      double lgx, lgy, lgxy;
      gsl_sf_lngamma_impl(x, &lgx);
      gsl_sf_lngamma_impl(y, &lgy);
      gsl_sf_lngamma_impl(x+y, &lgxy);
      *result = lgx + lgy - lgxy;
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_beta_impl(const double x, const double y, double * result)
{
  if(x < 50.0 && y < 50.0) {
    double gx, gy, gxy;
    gsl_sf_gamma_impl(x, &gx);
    gsl_sf_gamma_impl(y, &gy);
    gsl_sf_gamma_impl(x+y, &gxy);
    *result = (gx*gy)/gxy;
    return GSL_SUCCESS;
  }
  else {
    double lb;
    int stat_lb = gsl_sf_lnbeta_impl(x, y, &lb);
    if(stat_lb == GSL_SUCCESS) {
      return gsl_sf_exp_impl(lb, result);
    }
    else {
      *result = 0.0;
      return stat_lb;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_lnbeta_e(const double x, const double y, double * result)
{
  int status = gsl_sf_lnbeta_impl(x, y, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lnbeta_e", status);
  }
  return status;
}

int
gsl_sf_beta_e(const double x, const double y, double * result)
{
  int status = gsl_sf_beta_impl(x, y, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_beta_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_lnbeta(const double x, const double y)
{
  double r;
  int status = gsl_sf_lnbeta_impl(x, y, &r);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_lnbeta", status);
  }
  return r;
}

double
gsl_sf_beta(const double x, const double y)
{
  double r;
  int status = gsl_sf_beta_impl(x, y, &r);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_beta", status);
  }
  return r;
}
