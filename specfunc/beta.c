/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
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
      double lnpre_val;
      double lnpre_err;
      double lnpow_val;
      double lnpow_err;
      double t1, t2, t3;
      gsl_sf_result lnopr;
      gsl_sf_result gsx, gsy, gsxy;
      gsl_sf_gammastar_impl(x, &gsx);
      gsl_sf_gammastar_impl(y, &gsy);
      gsl_sf_gammastar_impl(x+y, &gsxy);
      gsl_sf_log_1plusx_impl(rat, &lnopr);
      lnpre_val = log(gsx.val*gsy.val/gsxy.val * M_SQRT2*M_SQRTPI);
      lnpre_err = gsx.err/gsx.val + gsy.err/gsy.val + gsxy.err/gsxy.val;
      t1 = min*log(rat);
      t2 = 0.5*log(min);
      t3 = (x+y-0.5)*lnopr.val;
      lnpow_val  = t1 - t2 - t3;
      lnpow_err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(t2) + fabs(t3));
      lnpow_err += fabs(x+y-0.5) * lnopr.err;
      result->val  = lnpre_val + lnpow_val;
      result->err  = lnpre_err + lnpow_err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result lgx, lgy, lgxy;
      int stat_gx  = gsl_sf_lngamma_impl(x, &lgx);
      int stat_gy  = gsl_sf_lngamma_impl(y, &lgy);
      int stat_gxy = gsl_sf_lngamma_impl(x+y, &lgxy);
      result->val  = lgx.val + lgy.val - lgxy.val;
      result->err  = lgx.err + lgy.err + lgxy.err;
      result->err += GSL_DBL_EPSILON * (fabs(lgx.val) + fabs(lgy.val) + fabs(lgxy.val));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_3(stat_gx, stat_gy, stat_gxy);
    }
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
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
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
