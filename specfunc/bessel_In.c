/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"

extern int gsl_sf_bessel_I1_scaled_impl(const double, double *);


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I_scaled_impl(int n, const double x, double * result)
{
  n = abs(n);  /* I(-n, z) = I(n, z) */
 
  if(n == 0) {
    *result = gsl_sf_bessel_I0_scaled(x);
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    double b;
    return gsl_sf_bessel_I1_scaled_impl(x, &b);
  }
  else { 
    if(x == 0.) {
      *result = 0.;
      return GSL_SUCCESS;
    }
    else {
      int j;
      int jmax = 2 * (n + (int) sqrt(50.*n));
      double renorm     = 1.e+12;
      double renorm_inv = 1.e-12;
      double two_over_x = 2./fabs(x);
      double b_np1 = 0.;
      double b_jp1 = 0.;
      double b_j   = 1.;
      double b_jm1;
    
      /* Downward recursion. [Gradshteyn + Ryzhik, 8.471.1] */
      for(j=jmax; j>0; j--){
      	b_jm1 = b_jp1 + j * two_over_x * b_j;
      	b_jp1 = b_j;
      	b_j   = b_jm1;
      
      	/* Renormalize to prevent overflow. */
      	if(fabs(b_j) > renorm){
	  b_np1 *= renorm_inv;
	  b_j   *= renorm_inv;
	  b_jp1 *= renorm_inv;
      	}
      
	/* Grab the value that we want on the way down. */
      	if(j == n) b_np1 = b_jp1;
      }

      /* Normalize using I0(x); j=0 here. */
      b_np1 *= gsl_sf_bessel_I0_scaled(x) / b_j;

      *result = x < 0. && (GSL_IS_ODD(n)) ? -b_np1 : b_np1;
      return GSL_SUCCESS;
    }
  }
}

int gsl_sf_bessel_I_impl(const int n, const double x, double * result)
{
  double ax = fabs(x);
  
  if(ax > GSL_LOG_DBL_MAX - 1.) {
    *result = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    double y;
    int status = gsl_sf_bessel_I_scaled_impl(n, x, &y);
    *result = exp(ax) * y;
    return status;
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_I_scaled_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_I_scaled_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_I_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_I_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_I_scaled(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_I_scaled_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I_scaled");
  }
  return y;
}

double gsl_sf_bessel_I(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_I_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I");
  }
  return y;
}
