/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"


/* i_{l+1}/i_l
 */
static
int
bessel_il_CF1(const int l, const double x, const double threshold, double * ratio)
{
  const int kmax = 2000;
  double tk   = 1.0;
  double sum  = 1.0;
  double rhok = 0.0;
  int k;

  for(k=1; k<=kmax; k++) {
    double ak = (x/(2.0*l+1.0+2.0*k)) * (x/(2.0*l+3.0+2.0*k));
    rhok = -ak*(1.0 + rhok)/(1.0 + ak*(1.0 + rhok));
    tk  *= rhok;
    sum += tk;
    if(fabs(tk/sum) < threshold) break;
  }

  *ratio = x/(2.0*l+3.0) * sum;

  if(k == kmax)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_i0_scaled_impl(const double x, gsl_sf_result * result)
{
  double ax = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(ax < 0.2) {
    const double eax = exp(-ax);
    const double y = ax*ax;
    const double c1 = 1.0/6.0;
    const double c2 = 1.0/120.0;
    const double c3 = 1.0/5040.0;
    const double c4 = 1.0/362880.0;
    const double c5 = 1.0/39916800.0;
    const double sum = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*c5))));
    result->val = eax * sum;
    result->err = 2.0 * GSL_DBL_EPSILON * result->val;
  }
  else if(ax < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = (1.0 - exp(-2.0*ax))/(2.0*ax);
    result->err = 2.0 * GSL_DBL_EPSILON * result->val;
  }
  else {
    result->val = 1.0/(2.0*ax);
    result->err = 2.0 * GSL_DBL_EPSILON * result->val;
  }
  return GSL_SUCCESS;
}


int gsl_sf_bessel_i1_scaled_impl(const double x, gsl_sf_result * result)
{
  double ax = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(ax < 3.0*GSL_DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(ax < 0.25) {
    const double eax = exp(-ax);
    const double y  = x*x;
    const double c1 = 1.0/10.0;
    const double c2 = 1.0/280.0;
    const double c3 = 1.0/15120.0;
    const double c4 = 1.0/1330560.0;
    const double c5 = 1.0/172972800.0;
    const double sum = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*c5))));
    result->val = eax * x/3.0 * sum;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double ex = exp(-2.0*ax);
    result->val = 0.5 * (ax*(1.0+ex) - (1.0-ex)) / (ax*ax);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    if(x < 0.0) result->val = -result->val;
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_i2_scaled_impl(const double x, gsl_sf_result * result)
{
  double ax = fabs(x);

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(ax < 4.0*GSL_SQRT_DBL_MIN) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(ax < 0.25) {
    const double y = x*x;
    const double c1 = 1.0/14.0;
    const double c2 = 1.0/504.0;
    const double c3 = 1.0/33264.0;
    const double c4 = 1.0/3459456.0;
    const double c5 = 1.0/518918400.0;
    const double sum = 1.0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*c5))));
    const double pre = exp(-ax) * x*x/15.0;
    result->val = pre * sum;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double ex = exp(-2.0*ax);
    double x2 = x*x;
    result->val = 0.5 * ((3.0+x2)*(1.0-ex) - 3.0*ax*(1.0+ex))/(ax*ax*ax);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}  


int gsl_sf_bessel_il_scaled_impl(const int l, double x, gsl_sf_result * result)
{
  double sgn = 1.0;
  double ax  = fabs(x);

  if(x < 0.0) {
    /* i_l(-x) = (-1)^l i_l(x) */
    sgn = ( GSL_IS_ODD(l) ? -1.0 : 1.0 );
    x = -x;
  }

  if(l < 0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 0) {
    gsl_sf_result il;
    int stat_il = gsl_sf_bessel_i0_scaled_impl(x, &il);
    result->val = sgn * il.val;
    result->err = il.err;
    return stat_il;
  }
  else if(l == 1) {
    gsl_sf_result il;
    int stat_il = gsl_sf_bessel_i1_scaled_impl(x, &il);
    result->val = sgn * il.val;
    result->err = il.err;
    return stat_il;
  }
  else if(l == 2) {
    gsl_sf_result il;
    int stat_il = gsl_sf_bessel_i2_scaled_impl(x, &il);
    result->val = sgn * il.val;
    result->err = il.err;
    return stat_il;
  }
  else if(x*x < 10.0*(l+1.5)/M_E) {
    gsl_sf_result b;
    int stat = gsl_sf_bessel_IJ_taylor_impl(l+0.5, x, 1, 50, GSL_DBL_EPSILON, &b);
    double pre   = exp(-ax) * sqrt((0.5*M_PI)/x);
    result->val  = sgn * pre * b.val;
    result->err  = pre * b.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat;
  }
  else if(l < 150) {
    gsl_sf_result i0_scaled;
    int stat_i0  = gsl_sf_bessel_i0_scaled_impl(ax, &i0_scaled);
    double rat;
    int stat_CF1 = bessel_il_CF1(l, ax, GSL_DBL_EPSILON, &rat);
    double iellp1 = rat * GSL_SQRT_DBL_MIN;
    double iell	  = GSL_SQRT_DBL_MIN;
    double iellm1;
    int ell;
    for(ell = l; ell >= 1; ell--) {
      iellm1 = iellp1 + (2*ell + 1)/x * iell;
      iellp1 = iell;
      iell   = iellm1;
    }
    result->val  = sgn * i0_scaled.val * (GSL_SQRT_DBL_MIN / iell);
    result->err  = i0_scaled.err * (GSL_SQRT_DBL_MIN / iell);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_i0, stat_CF1);
  }
  else if(GSL_MIN(0.29/(l*l+1.0), 0.5/(l*l+1.0+x*x)) < 0.5*GSL_ROOT3_DBL_EPSILON) {
    int status = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(l + 0.5, x, result);
    double pre = sqrt((0.5*M_PI)/x);
    result->val *= sgn * pre;
    result->err *= pre;
    return status;
  }
  else {
    /* recurse down from safe values */
    double rt_term = sqrt((0.5*M_PI)/x);
    const int LMAX = 2 + (int) (1.2 / GSL_ROOT6_DBL_EPSILON);
    gsl_sf_result r_iellp1;
    gsl_sf_result r_iell;
    int stat_a1 = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(LMAX + 1 + 0.5, x, &r_iellp1);
    int stat_a2 = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(LMAX     + 0.5, x, &r_iell);
    double iellp1 = r_iellp1.val;
    double iell   = r_iell.val;
    double iellm1;
    int ell;
    iellp1 *= rt_term;
    iell   *= rt_term;
    for(ell = LMAX; ell >= l+1; ell--) {
      iellm1 = iellp1 + (2*ell + 1)/x * iell;
      iellp1 = iell;
      iell   = iellm1;
    }
    result->val  = sgn * iellm1;
    result->err  = fabs(result->val)*(GSL_DBL_EPSILON + fabs(r_iellp1.err/r_iellp1.val) + fabs(r_iell.err/r_iell.val));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_a1, stat_a2);
  }
}


int gsl_sf_bessel_il_scaled_array_impl(const int lmax, const double x, double * result_array)
{
  int ell;
  gsl_sf_result r_iellp1;
  gsl_sf_result r_iell;
  int stat_0 = gsl_sf_bessel_il_scaled_impl(lmax+1, x, &r_iellp1);
  int stat_1 = gsl_sf_bessel_il_scaled_impl(lmax,   x, &r_iell);
  double iellp1 = r_iellp1.val;
  double iell   = r_iell.val;
  double iellm1;
  result_array[lmax] = iell;
  for(ell = lmax; ell >= 1; ell--) {
    iellm1 = iellp1 + (2*ell + 1)/x * iell;
    iellp1 = iell;
    iell   = iellm1;
    result_array[ell-1] = iellm1;
  }
  return GSL_ERROR_SELECT_2(stat_0, stat_1);
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_i0_scaled_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_i0_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_i0_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_i1_scaled_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_i1_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_i1_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_i2_scaled_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_i2_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_i2_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_il_scaled_e(const int l, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_il_scaled_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_il_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_il_scaled_array_e(const int lmax, const double x, double * il_array)
{
  int status = gsl_sf_bessel_il_scaled_array_impl(lmax, x, il_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_il_scaled_array_e", status);
  }
  return status;
}
