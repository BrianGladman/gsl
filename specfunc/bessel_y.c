/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* [Abramowitz+Stegun, 10.1.3]
 * with lmax=15, precision ~ 15D for x < 3
 *
 * checked OK [GJ] Wed May 13 15:41:25 MDT 1998 
 */
static int bessel_yl_small_x(int l, const double x, double * result)
{
  const int lmax = 15;
  int i;
  double num_fact;
  double den = gsl_sf_pow_int(x, l+1);
  if(gsl_sf_doublefact_impl(2*l-1, &num_fact) != GSL_SUCCESS || den == 0.0) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    double t = -0.5*x*x;
    double sum = 1.0;
    double t_coeff = 1.0;
    double t_power = 1.0;
    double delta;
    for(i=1; i<lmax; i++) {
      t_coeff /= i*(2*(i-l) - 1);
      t_power *= t;
      delta = t_power*t_coeff;
      sum += delta;
      if(fabs(delta) < GSL_MACH_EPS) break;
    }
    *result = -num_fact/den * sum;
    return GSL_SUCCESS;
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* checked OK [GJ] Wed May 13 15:19:15 MDT 1998 */
int gsl_sf_bessel_y0_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(1.0/DBL_MAX > 0.0 && x < 1.0/DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x < GSL_ROOT6_MACH_EPS) {
    double x2 = x*x;
    double x4 = x2*x2;
    *result = (-1.0 + 0.5*x2 - x4/24.)/x;
    return GSL_SUCCESS;
  }
  else {
    *result = -cos(x)/x;
    return GSL_SUCCESS;
  }
}

/* checked OK [GJ] Wed May 13 15:19:27 MDT 1998 */
int gsl_sf_bessel_y1_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(1.0/DBL_MAX > 0.0 && 144.0*x*x < 1.0/DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x < GSL_ROOT6_MACH_EPS) {
    double x2 = x*x;
    double x4 = x2*x2;
    *result = -(144.0 + 72.0*x2 - 18.0*x4)/(144.0*x2);
    return GSL_SUCCESS;
  }
  else {
    *result = -cos(x)/(x*x) - sin(x)/x;
    return GSL_SUCCESS;
  }
}

/* checked OK [GJ] Wed May 13 15:19:39 MDT 1998 */
int gsl_sf_bessel_y2_impl(const double x, double * result)
{
  if(x <= 0.0) {
    return GSL_EDOM;
  }
  else if(1.0/DBL_MAX > 0.0 && 1152.0*x*x*x < 1.0/DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(x < GSL_ROOT6_MACH_EPS) {
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x2*x2;
    *result = -(3456.0 + 576.0 * x2 + 144.0 * x4)/(1152.0 * x3);
    return GSL_SUCCESS;
  }
  else {
    double three_over_x2 = 3.0/(x*x);
    *result = (1.0 - three_over_x2)/x * cos(x) - three_over_x2 * sin(x);
    return GSL_SUCCESS;
  }
}

/* checked OK [GJ] Wed May 13 16:28:01 MDT 1998 */
int gsl_sf_bessel_yl_impl(int l, const double x, double * result)
{
  if(l < 0 || x <= 0.0) {
    return GSL_EDOM;
  }
  else if(x < 3.0) {
    return bessel_yl_small_x(l, x, result);
  }
  else if(GSL_ROOT3_MACH_EPS * x > (l*l + l + 1)) {
    int status = gsl_sf_bessel_Ynu_asympx_impl(l + 0.5, x, result);
    if(status == GSL_SUCCESS) *result *= sqrt(M_PI/(2.0*x));
    return status;
  }
  else if(l > 30) {
    int status = gsl_sf_bessel_Ynu_asymp_Olver_impl(l + 0.5, x, result);
    if(status == GSL_SUCCESS) *result *= sqrt(M_PI/(2.0*x));
    return status;
  }
  else if(l == 0) {
    return gsl_sf_bessel_y0_impl(x, result);
  }
  else if(l == 1) {
    return gsl_sf_bessel_y1_impl(x, result);
  }
  else if(l == 2) {
    return gsl_sf_bessel_y2_impl(x, result);
  }
  else {
    /* recurse upward */
    int j;
    double by, bym, byp;
    gsl_sf_bessel_y1_impl(x, &by);
    gsl_sf_bessel_y0_impl(x, &bym);
    for(j=1; j<l; j++) { 
      byp = (2*j+1)/x*by - bym;
      bym = by;
      by  = byp;
    }
    *result = by;
    return GSL_SUCCESS;
  }
}

/* checked OK [GJ] Wed May 13 16:33:10 MDT 1998 */
int gsl_sf_bessel_yl_array_impl(const int lmax, const double x, double * result_array)
{
  if(lmax < 1 || x <= 0.0) {
    return GSL_EDOM;
  }
  else {
    int ell;
    double yellp1, yell, yellm1;
    gsl_sf_bessel_y1_impl(x, &yell);
    gsl_sf_bessel_y0_impl(x, &yellm1);
    result_array[0] = yellm1;
    result_array[1] = yell;
    for(ell = 1; ell < lmax; ell++) {
      yellp1 = (2*ell+1)/x * yell - yellm1;
      result_array[ell+1] = yellp1;
      yellm1 = yell;
      yell   = yellp1;
    }
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_y0_e(const double x, double * result)
{
  int status = gsl_sf_bessel_y0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_y0_e", status);
  }
  return status;
}

int gsl_sf_bessel_y1_e(const double x, double * result)
{
  int status = gsl_sf_bessel_y1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_y1_e", status);
  }
  return status;
}

int gsl_sf_bessel_y2_e(const double x, double * result)
{
  int status = gsl_sf_bessel_y2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_y2_e", status);
  }
  return status;
}

int gsl_sf_bessel_yl_e(const int l, const double x, double * result)
{
  int status = gsl_sf_bessel_yl_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_yl_e", status);
  }
  return status;
}


int gsl_sf_bessel_yl_array_e(const int lmax, const double x, double * result_array)
{
  int status = gsl_sf_bessel_yl_array_impl(lmax, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_yl_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_y0(const double x)
{
  double y;
  int status = gsl_sf_bessel_y0_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_y0", status);
  }
  return y;
}

double gsl_sf_bessel_y1(const double x)
{
  double y;
  int status = gsl_sf_bessel_y1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_y1", status);
  }
  return y;
}

double gsl_sf_bessel_y2(const double x)
{
  double y;
  int status = gsl_sf_bessel_y2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_y2", status);
  }
  return y;
}

double gsl_sf_bessel_yl(const int l, const double x)
{
  double y;
  int status = gsl_sf_bessel_yl_impl(l, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_yl", status);
  }
  return y;
}
