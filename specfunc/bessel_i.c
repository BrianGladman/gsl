/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_i0_scaled_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < GSL_ROOT4_MACH_EPS) {
    *result = exp(-ax) * (1. + x*x/6.);
  }
  else {
    *result = (1.0 - exp(-2.0*ax))/(2.0*ax);
  }
  return GSL_SUCCESS;
}

int gsl_sf_bessel_i1_scaled_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < 3.*DBL_MIN) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else if(ax < 2.*GSL_ROOT4_MACH_EPS) {
    *result = exp(-ax) * x/3. * (1. + x*x/10.);
    return GSL_SUCCESS;
  }
  else {
    double ex = exp(-2.0*ax);
    *result = -0.5 * ((1.0-ex)/(ax*ax) - (1.0+ex)) / ax;
    return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_i2_scaled_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < GSL_SQRT_DBL_MIN) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else if(ax < 2.*GSL_ROOT4_MACH_EPS) {
    *result = exp(-ax) * x*x/15. * (1. + x*x/14.);
    return GSL_SUCCESS;
  }
  else {
    double ex = exp(-2.0*ax);
    double x2 = x*x;
    *result =  0.5 * ((3./x2 - 1.) * (1.0-ex)/ax - 3.*(1.0+ex)/x2);
    return GSL_SUCCESS;
  }
}  

int gsl_sf_bessel_il_scaled_impl(const int l, const double x, double * result)
{
  if(l < 0 || x < 0.0) {
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x*x < 10.*(l+1.5)*GSL_ROOT5_MACH_EPS) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_Jnu_taylor_impl(l+0.5, x, 1, 4, &b);
    *result = exp(-fabs(x)) * sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(GSL_ROOT3_MACH_EPS * x > (l*l + l + 1)) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_scaled_asympx_impl(l + 0.5, x, &b);
    *result = sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(l > 30) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(l + 0.5, x, &b);
    *result = sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(l == 0) {
    return gsl_sf_bessel_i0_scaled_impl(x, result);
  }
  else if(l == 1) {
    return gsl_sf_bessel_i1_scaled_impl(x, result);
  }
  else if(l == 2) {
    return gsl_sf_bessel_i2_scaled_impl(x, result);
  }
  else {
    /* recurse down from safe values */
    double rt_term = sqrt(M_PI/(2.0*x));
    double iellp1, iell, iellm1;
    const int LMAX = 31;
    int ell;
    gsl_sf_bessel_asymp_Inu_asymp_unif_impl(LMAX + 1 + 0.5, x, &iellp1);
    gsl_sf_bessel_asymp_Inu_asymp_unif_impl(LMAX     + 0.5, x, &iell);
    iellp1 *= rt_term;
    iell   *= rt_term;
    for(ell = LMAX; ell >= l+1; ell--) {
      iellm1 = iellp1 + (2*ell + 1)/x * iell;
      iellp1 = iell;
      iell   = iellm1;
    }
    *result = iellm1;
    return GSL_SUCCESS;
  }
}


int gsl_sf_bessel_il_scaled_array_impl(const int lmax, const double x, double * result_array)
{
  int ell;
  double iellp1, iell, iellm1;
  gsl_sf_bessel_il_scaled_impl(lmax+1, x, &iellp1);
  gsl_sf_bessel_il_scaled_impl(lmax,   x, &iell);
  result_array[lmax] = iell;
  for(ell = lmax; ell >= 1; ell--) {
    iellm1 = iellp1 + (2*ell + 1)/x * iell;
    iellp1 = iell;
    iell   = iellm1;
    result_array[ell-1] = iellm1;
  }
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_i0_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_i0_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_i0_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_i1_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_i1_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_i1_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_i2_scaled_e(const double x, double * result)
{
  int status = gsl_sf_bessel_i2_scaled_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_i2_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_il_scaled_e(const int l, const double x, double * result)
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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_i0_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_i0_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_i0_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_i1_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_i1_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_i1_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_i2_scaled(const double x)
{
  double y;
  int status = gsl_sf_bessel_i2_scaled_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_i2_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_il_scaled(const int l, const double x)
{
  double y;
  int status = gsl_sf_bessel_il_scaled_impl(l, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_il_scaled", status);
  }
  return y;
}
