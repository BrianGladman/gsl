/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"

#include "gamma_impl.h"
#include "bessel.h"
#include "bessel_K0_impl.h"
#include "bessel_K1_impl.h"

#include "bessel_Kn_impl.h"

#define Min(a,b) ((a) < (b) ? (a) : (b))

/* assumes: n >= 2
 * choose x*x < 4 MACH_EPS^(1/2) / n
 */
static int bessel_Kn_scaled_taylor(const int n, const double x, double * result)
{
  double ln_x_2 = log(0.5*x);
  double lngamn;
  double ln_pre;
  gsl_sf_lngamma_impl(n, &lngamn);
  ln_pre = lngamn - n*ln_x_2;
  if(ln_pre > GSL_LOG_DBL_MAX) {
    *result = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(ln_pre < GSL_LOG_DBL_MIN + 1.) {
    *result = 0.;
    return GSL_EUNDRFLW;
  }
  else {
    double ep = exp(ln_pre);
    double sgn = (GSL_IS_ODD(n) ? 1. : -1. );
    double sum = 1. - 0.25*x*x/(n-1.);
    double term1 = 0.5 * ep * sum;
    double term2 = sgn/(n*ep) * ln_x_2;
    *result = exp(x) * (term1 + term2);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Kn_scaled_impl(int n, const double x, double * result)
{
  n = abs(n); /* K(-n, z) = K(n, z) */
  
  if(x <= 0.0) return GSL_EDOM;

  if(n == 0) {
    return gsl_sf_bessel_K0_scaled_impl(x, result);
  }
  else if(n == 1) {
    return gsl_sf_bessel_K1_scaled_impl(x, result);
  }
  else if(x*x * n < 4. * GSL_SQRT_MACH_EPS) {
    return bessel_Kn_scaled_taylor(n, x, result);
  }
  else if(GSL_ROOT3_MACH_EPS * x > 0.25 * (n*n + 1)) {
    return gsl_sf_bessel_Knu_scaled_asympx_impl((double)n, x, result);
  }
  else if(Min( 0.29/(n*n), 0.5/(n*n + x*x) ) < GSL_ROOT3_MACH_EPS) {
    return gsl_sf_bessel_Knu_scaled_asymp_unif_impl((double)n, x, result);
  }
  else {
    /* FIXME: must check for overflow */
    /* Upward recurrence. [Gradshteyn + Ryzhik, 8.471.1] */
    int j;
    double two_over_x = 2.0/x;
    double b_jm1;
    double b_j;
    double b_jp1;
    gsl_sf_bessel_K0_scaled_impl(x, &b_jm1);
    gsl_sf_bessel_K1_scaled_impl(x, &b_j);

    for(j=1; j<n; j++) {
      b_jp1 = b_jm1 + j * two_over_x * b_j;
      b_jm1 = b_j;
      b_j   = b_jp1; 
    } 
    
    *result = b_j; 
    return GSL_SUCCESS;
  }
}

int gsl_sf_bessel_Kn_scaled_array_impl(int n, const double x, double * result_array)
{
  n = abs(n); /* K(-n, z) = K(n, z) */
  
  if(n == 0) {
    return gsl_sf_bessel_K0_scaled_impl(x, &(result_array[0]));
  }
  else if(n == 1) {
    gsl_sf_bessel_K0_scaled_impl(x, &(result_array[0]));
    return gsl_sf_bessel_K1_scaled_impl(x, &(result_array[1]));
  }
  else {
    if(x <= 0.) {
      return GSL_EDOM;
    }
    else {
      /* Upward recurrence. [Gradshteyn + Ryzhik, 8.471.1] */
      int j;
      double two_over_x = 2.0/x;
      double b_jm1;
      double b_j;
      double b_jp1;
      gsl_sf_bessel_K0_scaled_impl(x, &b_jm1);
      gsl_sf_bessel_K1_scaled_impl(x, &b_j);
      result_array[0] = b_jm1;

      for(j=1; j<n; j++) {
	b_jp1 = b_jm1 + j * two_over_x * b_j;
	
	result_array[j]   = b_j;
	result_array[j+1] = b_jp1;

	b_jm1 = b_j;
	b_j   = b_jp1; 
      }

      return GSL_SUCCESS;
    }
  }
}

void testy(void)
{
  int i;
  double x = 2.e-4;
  double y;
  
  for(i=0; i<70; i++) {
    int stat = gsl_sf_bessel_Kn_scaled_impl(i, x, &y);
    printf("%3d  %26.16g  %26.16g   %d\n", i, x, y, stat);
  }
  exit(0);
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Kn_scaled_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_Kn_scaled_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Kn_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_Kn_e(const int n, const double x, double * result)
{
  double y = 0.;
  int status = gsl_sf_bessel_Kn_scaled_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Kn_e", status);
  }
  *result = exp(-x) * y;
  return status;
}

int gsl_sf_bessel_Kn_scaled_array_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_Kn_scaled_array_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Kn_scaled_array_e", status);
  }
  return status;
}

int gsl_sf_bessel_Kn_array_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_Kn_scaled_array_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Kn_array_e", status);
  }
  else {
    int i;
    double ex = exp(-x);
    for(i=0; i<=n; i++) result[i] *= ex;
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_Kn_scaled(const int n, const double x)
{
  double y = 0.;
  int status = gsl_sf_bessel_Kn_scaled_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Kn_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_Kn(const int n, const double x)
{
  double y = 0.;
  int status = gsl_sf_bessel_Kn_scaled_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Kn", status);
  }
  y *= exp(-x);
  return y;
}
