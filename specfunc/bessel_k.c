/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_bessel.h"



/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* [Abramowitz+Stegun, 10.2.4 + 10.2.6]
 * with lmax=15, precision ~ 15D for x < 3
 *
 * assumes l >= 1
 */
static int bessel_kl_scaled_small_x(int l, const double x, gsl_sf_result * result)
{
  gsl_sf_result num_fact;
  double den  = gsl_sf_pow_int(x, l+1);
  int stat_df = gsl_sf_doublefact_impl((unsigned int) (2*l-1), &num_fact);

  if(stat_df != GSL_SUCCESS || den == 0.0) {
    result->val = 0.0; /* FIXME: should be Inf */
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    const int lmax = 50;
    gsl_sf_result ipos_term;
    double ineg_term;
    double sgn = (GSL_IS_ODD(l) ? -1.0 : 1.0);
    double ex  = exp(x);
    double t = 0.5*x*x;
    double sum = 1.0;
    double t_coeff = 1.0;
    double t_power = 1.0;
    double delta;
    int stat_il;
    int i;

    for(i=1; i<lmax; i++) {
      t_coeff /= i*(2*(i-l) - 1);
      t_power *= t;
      delta = t_power*t_coeff;
      sum += delta;
      if(fabs(delta/sum) < GSL_DBL_EPSILON) break;
    }

    stat_il = gsl_sf_bessel_il_scaled_impl(l, x, &ipos_term);
    ineg_term =  sgn * num_fact.val/den * sum;
    result->val = -sgn * 0.5*M_PI * (ex*ipos_term.val - ineg_term);
    result->val *= ex;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_il;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_k0_scaled_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    result->val = M_PI/(2.0*x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    if(result->val == 0.0)
      return GSL_EUNDRFLW;
    else 
      return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_k1_scaled_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x < (M_SQRTPI+1.0)/(M_SQRT2*GSL_SQRT_DBL_MAX)) {
    result->val = 0.0;  /* FIXME: should be Inf */
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    result->val = M_PI/(2.0*x) * (1.0 + 1.0/x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    if(result->val == 0.0)
      return GSL_EUNDRFLW;
    else 
      return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_k2_scaled_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(x < 2.0/GSL_ROOT3_DBL_MAX) {
    result->val = 0.0; /* FIXME: should be Inf */
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    result->val = M_PI/(2.0*x) * (1.0 + 3.0/x * (1.0 + 1.0/x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    if(result->val == 0.0)
      return GSL_EUNDRFLW;
    else 
      return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_kl_scaled_impl(int l, const double x, gsl_sf_result * result)
{
  if(l < 0 || x <= 0.0) {
    return GSL_EDOM;
  }
  else if(l == 0) {
    return gsl_sf_bessel_k0_scaled_impl(x, result);
  }
  else if(l == 1) {
    return gsl_sf_bessel_k1_scaled_impl(x, result);
  }
  else if(l == 2) {
    return gsl_sf_bessel_k2_scaled_impl(x, result);
  }
  else if(x < 3.0) {
    return bessel_kl_scaled_small_x(l, x, result);
  }
  else if(GSL_ROOT3_DBL_EPSILON * x > (l*l + l + 1)) {
    int status = gsl_sf_bessel_Knu_scaled_asympx_impl(l + 0.5, x, result);
    double pre = sqrt((0.5*M_PI)/x);
    result->val *= pre;
    result->err *= pre;
    return status;
  }
  else if(GSL_MIN(0.29/(l*l+1.0), 0.5/(l*l+1.0+x*x)) < GSL_ROOT3_DBL_EPSILON) {
    int status = gsl_sf_bessel_Knu_scaled_asymp_unif_impl(l + 0.5, x, result);
    double pre = sqrt((0.5*M_PI)/x);
    result->val *= pre;
    result->err *= pre;
    return status;
  }
  else {
    /* recurse upward */
    gsl_sf_result r_bk;
    gsl_sf_result r_bkm;
    int stat_1 = gsl_sf_bessel_k1_scaled_impl(x, &r_bk);
    int stat_0 = gsl_sf_bessel_k0_scaled_impl(x, &r_bkm);
    double bkp;
    double bk  = r_bk.val;
    double bkm = r_bkm.val;
    int j;
    for(j=1; j<l; j++) { 
      bkp = (2*j+1)/x*bk + bkm;
      bkm = bk;
      bk  = bkp;
    }
    result->val  = bk;
    result->err  = fabs(bk) * (fabs(r_bk.err/r_bk.val) + fabs(r_bkm.err/r_bkm.val));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_1, stat_0);
  }
}

int gsl_sf_bessel_kl_scaled_array_impl(const int lmax, const double x, double * result_array)
{
  if(lmax < 1 || x <= 0.0) {
    return GSL_EDOM;
  }
  else {
    int ell;
    double kellp1, kell, kellm1;
    gsl_sf_result r_kell;
    gsl_sf_result r_kellm1;
    gsl_sf_bessel_k1_scaled_impl(x, &r_kell);
    gsl_sf_bessel_k0_scaled_impl(x, &r_kellm1);
    kell   = r_kell.val;
    kellm1 = r_kellm1.val;
    result_array[0] = kellm1;
    result_array[1] = kell;
    for(ell = 1; ell < lmax; ell++) {
      kellp1 = (2*ell+1)/x * kell + kellm1;
      result_array[ell+1] = kellp1;
      kellm1 = kell;
      kell   = kellp1;
    }
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_k0_scaled_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_k0_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_k0_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_k1_scaled_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_k1_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_k1_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_k2_scaled_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_k2_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_k2_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_kl_scaled_e(const int l, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_kl_scaled_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_kl_scaled_e", status);
  }
  return status;
}


int gsl_sf_bessel_kl_scaled_array_e(const int lmax, const double x, double * result_array)
{
  int status = gsl_sf_bessel_kl_scaled_array_impl(lmax, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_kl_scaled_array_e", status);
  }
  return status;
}
