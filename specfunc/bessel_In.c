/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"



extern int gsl_sf_bessel_I1_scaled_impl(const double, double *);



/* n >= 2; x >= 0
 * checked OK [GJ]
 */
static int bessel_In_scaled(const int n, const double x, double * b_n, double * b_nm1)
{
  if(x == 0.) {
    *b_n   = 0.;
    *b_nm1 = 0.;
    return GSL_SUCCESS;
  }
  else if(x*x < 4.*(n+1)*GSL_SQRT_MACH_EPS) {
    double ex = exp(-x);
    gsl_sf_bessel_Inu_Jnu_taylor_impl(n-1, x, 1, 4, b_nm1);
    gsl_sf_bessel_Inu_Jnu_taylor_impl(n,   x, 1, 4, b_n  );
    *b_nm1 *= ex;
    *b_n   *= ex;
    if(*b_n == 0. || *b_nm1 == 0.) {
      return GSL_EUNDRFLW;
    }
    else {
      return GSL_SUCCESS;
    }
  }
  else {
    int j;
    int jmax = 2 * (n + (int) sqrt(50.*n));
    double renorm     = 1.e+12;
    double renorm_inv = 1.e-12;
    double two_over_x = 2./fabs(x);
    double ratio;
    double b_jp1 = 0.;
    double b_j   = 1.;
    double b_jm1;
    *b_n   = 0.;
    *b_nm1 = 0.;
  
    /* backward recursion [Gradshteyn + Ryzhik, 8.471.1] */
    for(j=jmax; j>0; j--){
      b_jm1 = b_jp1 + j * two_over_x * b_j;
    
      /* Renormalize to prevent overflow. */
      if(fabs(b_jm1) > renorm){
        *b_nm1 *= renorm_inv;
        *b_n   *= renorm_inv;
        b_j   *= renorm_inv;
        b_jm1 *= renorm_inv;
      }
    
      /* Grab the values that we want on the way down. */
      if(j == n) {
        *b_n   = b_j;
	*b_nm1 = b_jm1;
      }

      b_jp1 = b_j;
      b_j   = b_jm1;
    }

    /* Normalize using I0(x); j=0 here. */
    ratio = gsl_sf_bessel_I0_scaled(x) / b_j;
    *b_nm1 *= ratio;
    *b_n   *= ratio;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int gsl_sf_bessel_In_scaled_impl(int n, const double x, double * result, double * harvest)
{
  n = abs(n);  /* I(-n, z) = I(n, z) */
 
  if(n == 0) {
    *result = gsl_sf_bessel_I0_scaled(x);
    if(harvest != 0) harvest[0] = *result;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    double b;
    int status = gsl_sf_bessel_I1_scaled_impl(x, &b);
    if(harvest != 0) {
      harvest[0] = gsl_sf_bessel_I0_scaled(x);
      harvest[1] = b;
    }
    *result = b;
    return GSL_SUCCESS;
  }
  else { 
    if(x == 0.) {
      int j;
      *result = 0.;
      if(harvest != 0) {
        harvest[0] = 1.;
        for(j=1; j<=n; j++) harvest[j] = 0.;
      }
      return GSL_SUCCESS;
    }
    else {
      int j;
      int jmax = 2 * (n + (int) sqrt(50.*n));
      double renorm     = 1.e+12;
      double renorm_inv = 1.e-12;
      double two_over_x = 2./fabs(x);
      double ratio;
      double b_n;
      double b_np1 = 0.;
      double b_jp1 = 0.;
      double b_j   = 1.;
      double b_jm1;
    
      /* backward recursion [Gradshteyn + Ryzhik, 8.471.1] */
      for(j=jmax; j>0; j--){
      	b_jm1 = b_jp1 + j * two_over_x * b_j;
      	b_jp1 = b_j;
      	b_j   = b_jm1;
      
      	/* Renormalize to prevent overflow. */
      	if(fabs(b_j) > renorm){
	  b_np1 *= renorm_inv;
	  b_n   *= renorm_inv;
	  b_j   *= renorm_inv;
	  b_jp1 *= renorm_inv;
      	}
      
	/* Grab the value(s) that we want on the way down. */
      	if(j == n)   b_np1 = b_jp1;
	if(j == n-1) b_n   = b_jp1;
      }

      /* Normalize using I0(x); j=0 here. */
      ratio = gsl_sf_bessel_I0_scaled(x) / b_j;
      b_np1 *= ratio;
      b_n   *= ratio;
      *result = ( x < 0. && (GSL_IS_ODD(n)) ? -b_np1 : b_np1 );
      
      if(harvest != 0) {
        harvest[n]   = *result;
	harvest[n-1] = ( x < 0. && (GSL_IS_ODD(n-1)) ? -b_n : b_n );
        b_jp1 = harvest[n];
	b_j   = harvest[n-1];
        for(j=n-2; j>=0; j--){
      	  b_jm1 = b_jp1 + j * two_over_x * b_j;
      	  b_jp1 = b_j;
      	  b_j   = b_jm1;
	  harvest[j] = ( x < 0. && (GSL_IS_ODD(j)) ? b_jm1 : -b_jm1 );
        }
      }

      return GSL_SUCCESS;
    }
  }
}

int gsl_sf_bessel_In_impl(const int n, const double x, double * result, double * harvest)
{
  int j;
  double ax = fabs(x);
  
  if(ax > GSL_LOG_DBL_MAX - 1.) {
    *result = 0.; /* FIXME: should be Inf */
    if(harvest != 0) {
      for(j=0; j<=n; j++) harvest[j] = 0.; /* FIXME: should be Inf */
    }
    return GSL_EOVRFLW;
  }
  else {
    double y;
    double eax = exp(ax);
    int status = gsl_sf_bessel_In_scaled_impl(n, x, &y, harvest);
    *result = eax * y;
    if(harvest != 0) {
      for(j=0; j<=n; j++) harvest[j] *= eax;
    }
    return status;
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_In_scaled_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_In_scaled_impl(n, x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_In_impl(n, x, result, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_scaled_array_e(const int n, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_bessel_In_scaled_impl(n, x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I_scaled_array_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_array_e(const int n, const double x, double * result_array)
{
  double y;
  int status = gsl_sf_bessel_In_impl(n, x, &y, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_I_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_In_scaled(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_In_scaled_impl(n, x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I_scaled");
  }
  return y;
}

double gsl_sf_bessel_In(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_In_impl(n, x, &y, (double *)0);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_I");
  }
  return y;
}
