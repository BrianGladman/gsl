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
static int bessel_kl_scaled_small_x(int l, const double x, double * result)
{
  const int lmax = 15;
  int i;
  double num_fact;
  double den = gsl_sf_pow_int(x, l+1);
  if(gsl_sf_doublefact_impl((unsigned int) 2*l-1, &num_fact) != GSL_SUCCESS || den == 0.0) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    double ipos_term;
    double ineg_term;
    double sgn = (GSL_IS_ODD(l) ? -1.0 : 1.0);
    double ex  = exp(x);
    double t = 0.5*x*x;
    double sum = 1.0;
    double t_coeff = 1.0;
    double t_power = 1.0;
    double delta;
    for(i=1; i<lmax; i++) {
      t_coeff /= i*(2*(i-l) - 1);
      t_power *= t;
      delta = t_power*t_coeff;
      sum += delta;
      if(fabs(delta/sum) < GSL_MACH_EPS) break;
    }
    gsl_sf_bessel_il_scaled_impl(l, x, &ipos_term);
    ineg_term =  sgn * num_fact/den * sum;
    *result   = -sgn * 0.5*M_PI * (ex*ipos_term - ineg_term);
    *result  *= ex;
    return GSL_SUCCESS;
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_k0_scaled_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else {
    *result = M_PI/(2.0*x);
    if(*result == 0.0)
      return GSL_EUNDRFLW;
    else 
      return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_k1_scaled_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(1.0/DBL_MAX > 0.0 && 2.0*x*x < M_PI/DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    *result = M_PI/(2.0*x) * (1.0 + 1.0/x);
    if(*result == 0.0)
      return GSL_EUNDRFLW;
    else 
      return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_k2_scaled_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(1.0/DBL_MAX > 0.0 && 2.0*x*x*x < 3.0*M_PI/DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    *result = M_PI/(2.0*x) * (1.0 + 3.0/x * (1.0 + 1.0/x));
    if(*result == 0.0)
      return GSL_EUNDRFLW;
    else 
      return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_kl_scaled_impl(int l, const double x, double * result)
{
  if(l < 0 || x <= 0.0) {
    return GSL_EDOM;
  }
  else if(x < 3.0) {
    return bessel_kl_scaled_small_x(l, x, result);
  }
  else if(GSL_ROOT3_MACH_EPS * x > (l*l + l + 1)) {
    int status = gsl_sf_bessel_Knu_scaled_asympx_impl(l + 0.5, x, result);
    if(status == GSL_SUCCESS) *result *= sqrt(M_PI/(2.0*x));
    return status;
  }
  else if(GSL_MIN(0.29/(l*l+1.), 0.5/(l*l+1.+x*x)) < GSL_ROOT3_MACH_EPS) {
    int status = gsl_sf_bessel_Knu_scaled_asymp_unif_impl(l + 0.5, x, result);
    if(status == GSL_SUCCESS) *result *= sqrt(M_PI/(2.0*x));
    return status;
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
  else {
    /* recurse upward */
    int j;
    double bk, bkm, bkp;
    gsl_sf_bessel_k1_scaled_impl(x, &bk);
    gsl_sf_bessel_k0_scaled_impl(x, &bkm);
    for(j=1; j<l; j++) { 
      bkp = (2*j+1)/x*bk + bkm;
      bkm = bk;
      bk  = bkp;
    }
    *result = bk;
    return GSL_SUCCESS;
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
    gsl_sf_bessel_k1_scaled_impl(x, &kell);
    gsl_sf_bessel_k0_scaled_impl(x, &kellm1);
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

int gsl_sf_bessel_k0_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_k0_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_k0_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_k1_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_k1_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_k1_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_k2_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_k2_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_k2_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_kl_scaled_e(const int l, const double x, double * result)
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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_k0_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_k0_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_k0_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_k1_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_k1_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_k1_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_k2_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_k2_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_k2_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_kl_scaled(const int l, const double x)
{
  double y;
  int status = gsl_sf_bessel_kl_scaled_impl(l, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_kl_scaled", status);
  }
  return y;
}
