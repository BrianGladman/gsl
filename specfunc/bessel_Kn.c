/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"

extern int gsl_sf_bessel_K0_scaled_impl(const double, double *);
extern int gsl_sf_bessel_K1_scaled_impl(const double, double *);


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
  else {
    if(x <= 0.) {
      return GSL_EDOM;
    }
    else {
      /* Upward recurrence. [Gradshteyn + Ryzhik, 8.471.1] */
      int j;
      double two_over_x = 2.0/x;
      double b_jm1 = gsl_sf_bessel_K0_scaled(x);
      double b_j   = gsl_sf_bessel_K1_scaled(x);
      double b_jp1;

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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_K_scaled_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_K_scaled_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_K_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_K_e(const int n, const double x, double * result)
{
  double y = 0.;
  int status = gsl_sf_bessel_K_scaled_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_K_e", status);
  }
  *result = exp(-x) * y;
  return status;
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_K(const int n, const double x)
{
  double y = 0.;
  int status = gsl_sf_bessel_K_scaled_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_K");
  }
  y *= exp(-x);
  return y;
}
