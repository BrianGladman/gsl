/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"

#include "bessel.h"
#include "bessel_K0_impl.h"
#include "bessel_K1_impl.h"

#include "bessel_Kn_impl.h"


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Kn_scaled_impl(int n, const double x, double * result)
{
  n = abs(n); /* K(-n, z) = K(n, z) */
  
  if(n == 0) {
    return gsl_sf_bessel_K0_scaled_impl(x, result);
  }
  else if(n == 1) {
    return gsl_sf_bessel_K1_scaled_impl(x, result);
  }
  else if(x*x < 4.*(n+1)*GSL_SQRT_MACH_EPS) {
    return gsl_sf_bessel_Knu_asympx_impl((double)n, x, result);
  }
  else if(n > 700) {
    return gsl_sf_bessel_Knu_scaled_asymp_unif_impl((double)n, x, result);
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

      for(j=1; j<n; j++) {
	b_jp1 = b_jm1 + j * two_over_x * b_j;
	b_jm1 = b_j;
	b_j   = b_jp1; 
      } 
      
      *result = b_j; 
      return GSL_SUCCESS;
    }
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
  int status = gsl_sf_bessel_K_scaled_impl(n, x, &y);
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
  int i;
  double ex = exp(-x);
  int status = gsl_sf_bessel_Kn_scaled_array_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Kn_array_e", status);
  }
  else {
    for(i=0; i<=n; i++) result[i] *= exp(-x);
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
