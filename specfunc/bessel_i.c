/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"

#define locMin(a,b) ((a) < (b) ? (a) : (b))
#define locMax(a,b) ((a) > (b) ? (a) : (b))


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_i0_scaled_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < GSL_ROOT4_MACH_EPS) {
    *result = exp(-ax) * (1.0 + x*x/6.0);
  }
  else if(ax < -0.5*GSL_LOG_MACH_EPS){
    *result = (1.0 - exp(-2.0*ax))/(2.0*ax);
  }
  else {
    *result = 1.0/(2.0*ax);
  }
  return GSL_SUCCESS;
}

int gsl_sf_bessel_i1_scaled_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < 3.0*DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(ax < 0.02) {
    *result = exp(-ax) * x/3.0 * (1.0 + x*x/10.0 * (1.0 + x*x/28.0 * (1.0 + x*x/54.0)));
    return GSL_SUCCESS;
  }
  else {
    double ex = exp(-2.0*ax);
    *result = 0.5 * (ax*(1.0+ex) - (1.0-ex)) / (ax*ax);
    if(x < 0.0) *result = - *result;
    return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_i2_scaled_impl(const double x, double * result)
{
  double ax = fabs(x);

  if(ax < 4.0*GSL_SQRT_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else if(ax < 0.25) {
    double pre = exp(-ax) * x*x/15.0;
    *result = pre * (1.0 + x*x/14.0 * (1.0 + x*x/36.0 * (1.0 + x*x/66.0 * (1.0 + x*x/104.0 * (1.0 + x*x/150.0)))));
    return GSL_SUCCESS;
  }
  else {
    double ex = exp(-2.0*ax);
    double x2 = x*x;
    *result =  0.5 * ((3.0+x2)*(1.0-ex) - 3.0*ax*(1.0+ex))/(ax*ax*ax);
    return GSL_SUCCESS;
  }
}  

int gsl_sf_bessel_il_scaled_impl(const int l, double x, double * result)
{
  double sgn = 1.0;
  double ax  = fabs(x);

  if(x < 0.0) {
    /* i_l(-x) = (-1)^l i_l(x) */
    sgn = ( GSL_IS_ODD(l) ? -1.0 : 1.0 );
    x = -x;
  }

  if(l < 0) {
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(x*x < 10.0*(l+1.5)*GSL_ROOT5_MACH_EPS) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_Jnu_taylor_impl(l+0.5, x, 1, 4, &b);
    *result = sgn * exp(-ax) * sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(x*x < 10.0*(l+1.5)/M_E) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_Jnu_taylor_impl(l+0.5, x, 1, 14, &b);
    *result = sgn * exp(-ax) * sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(GSL_ROOT3_MACH_EPS * x > (l*l + l + 1)) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_scaled_asympx_impl(l + 0.5, x, &b);
    *result = sgn * sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(locMin(0.29/(l*l+1.), 0.5/(l*l+1.+x*x)) < GSL_ROOT3_MACH_EPS) {
    double b = 0.0;
    int status = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(l + 0.5, x, &b);
    *result = sgn * sqrt(M_PI/(2.0*x)) * b;
    return status;
  }
  else if(l == 0) {
    double il;
    int stat_il = gsl_sf_bessel_i0_scaled_impl(x, &il);
    *result = sgn * il;
    return stat_il;
  }
  else if(l == 1) {
    double il;
    int stat_il = gsl_sf_bessel_i1_scaled_impl(x, &il);
    *result = sgn * il;
    return stat_il;
  }
  else if(l == 2) {
    double il;
    int stat_il = gsl_sf_bessel_i2_scaled_impl(x, &il);
    *result = sgn * il;
    return stat_il;
  }
  else {
    /* recurse down from safe values */
    double rt_term = sqrt(M_PI/(2.0*x));
    double iellp1, iell, iellm1;
    int ell;
    const int LMAX = sqrt(locMax(0.5/GSL_ROOT3_MACH_EPS - x*x,
                          0.29/GSL_ROOT3_MACH_EPS));
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl(LMAX + 1 + 0.5, x, &iellp1);
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl(LMAX     + 0.5, x, &iell);
    iellp1 *= rt_term;
    iell   *= rt_term;
    for(ell = LMAX; ell >= l+1; ell--) {
      iellm1 = iellp1 + (2*ell + 1)/x * iell;
      iellp1 = iell;
      iell   = iellm1;
    }
    *result = sgn * iellm1;
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
