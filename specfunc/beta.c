/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_log.h"
#include "gsl_sf_psi.h"
#include "gsl_sf_gamma.h"


int
gsl_sf_lnbeta_impl(const double x, const double y, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0 || y <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    const double max = GSL_MAX(x,y);
    const double min = GSL_MIN(x,y);
    const double rat = min/max;

    if(rat < 0.2) {
      /* min << max, so be careful
       * with the subtraction
       */
      double lnpow;
      double lnpre;
      gsl_sf_result lnopr;
      gsl_sf_result gsx, gsy, gsxy;
      gsl_sf_gammastar_impl(x, &gsx);
      gsl_sf_gammastar_impl(y, &gsy);
      gsl_sf_gammastar_impl(x+y, &gsxy);
      gsl_sf_log_1plusx_impl(rat, &lnopr);
      lnpre = log(gsx.val*gsy.val/gsxy.val * M_SQRT2*M_SQRTPI);
      lnpow = min*log(rat) - 0.5*log(min) - (x+y-0.5)*lnopr.val;
      result->val = lnpre + lnpow;
      result->err = GSL_DBL_EPSILON * (fabs(lnpre) + fabs(lnpow));
    }
    else {
      gsl_sf_result lgx, lgy, lgxy;
      gsl_sf_lngamma_impl(x, &lgx);
      gsl_sf_lngamma_impl(y, &lgy);
      gsl_sf_lngamma_impl(x+y, &lgxy);
      result->val = lgx.val + lgy.val - lgxy.val;
      result->err = lgx.err + lgy.err + lgxy.err;
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_beta_impl(const double x, const double y, gsl_sf_result * result)
{
  if(x < 50.0 && y < 50.0) {
    gsl_sf_result gx, gy, gxy;
    gsl_sf_gamma_impl(x, &gx);
    gsl_sf_gamma_impl(y, &gy);
    gsl_sf_gamma_impl(x+y, &gxy);
    result->val  = (gx.val*gy.val)/gxy.val;
    result->err  = gx.err * gy.val/gxy.val;
    result->err += gy.err * gx.val/gxy.val;
    result->err += (gx.val*gy.val)/(gxy.val*gxy.val) * gxy.err;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result lb;
    int stat_lb = gsl_sf_lnbeta_impl(x, y, &lb);
    if(stat_lb == GSL_SUCCESS) {
      return gsl_sf_exp_err_impl(lb.val, lb.err, result);
    }
    else {
      result->val = 0.0;
      result->err = 0.0;
      return stat_lb;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_lnbeta_e(const double x, const double y, gsl_sf_result * result)
{
  int status = gsl_sf_lnbeta_impl(x, y, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_lnbeta_e", status);
  }
  return status;
}

int
gsl_sf_beta_e(const double x, const double y, gsl_sf_result * result)
{
  int status = gsl_sf_beta_impl(x, y, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_beta_e", status);
  }
  return status;
}
