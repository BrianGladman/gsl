/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"

#include "bessel.h"

#define Max(a,b) ((a) > (b) ? (a) : (b))
#define Min(a,b) ((a) < (b) ? (a) : (b))


/* standard backward recursion with arbitrary starting values
 * assumes: n >= 2, x > 0
 * checked OK [GJ] Sun May  3 09:09:40 EDT 1998 
 */
static void recurse_In(const int n, const double x, double I0, double * b_n, double * b_nm1)
{
  int j;
  const int jmax = 2 * (n + (int) sqrt(50.*n));
  const double renorm	  = 1.e+8;
  const double renorm_inv = 1.e-8;
  const double two_over_x = 2./x;
  double ratio;
  double b_jp1 = 0.;
  double b_j   = 1.;
  double b_jm1;
  double local_b_n   = 0.;
  double local_b_nm1 = 0.;

  /* backward recursion [Gradshteyn + Ryzhik, 8.471.1] */
  for(j=jmax; j>0; j--){
    b_jm1 = b_jp1 + j * two_over_x * b_j;
  
    /* Renormalize to prevent overflow. */
    if(fabs(b_jm1) > renorm){
      local_b_nm1 *= renorm_inv;
      local_b_n   *= renorm_inv;
      b_j   *= renorm_inv;
      b_jm1 *= renorm_inv;
    }
  
    /* Grab the values that we want on the way down. */
    if(j == n) {
      local_b_n   = b_j;
      local_b_nm1 = b_jm1;
    }

    b_jp1 = b_j;
    b_j   = b_jm1;
  }

  /* Normalize using I0(x); j=0 here. */
  ratio  = I0 / b_j;
  *b_n   = local_b_n   * ratio;
  *b_nm1 = local_b_nm1 * ratio;
}


/* backward recursion using asymptotic starting values
 * assumes: n >= 2, x > 0
 * checked OK [GJ] Sun May  3 09:15:06 EDT 1998 
 */
static void asymp_recurse_In_scaled(const int n, const double x, double * b_n, double * b_nm1)
{
  int j;
  double b_jp1, b_j, b_jm1;
  double two_over_x = 2./x;
  
  /* go high enough to apply asymptotic result */
  int big_nu = 2 + Max(n, (int)sqrt(Max(0., 0.5/GSL_ROOT3_MACH_EPS - x*x)));

  gsl_sf_bessel_Inu_scaled_asymp_unif_impl(big_nu  , x, &b_jp1);
  gsl_sf_bessel_Inu_scaled_asymp_unif_impl(big_nu-1, x, &b_j);

  for(j=big_nu-1; j>=n; j--){
    b_jm1 = b_jp1 + j * two_over_x * b_j;
    b_jp1 = b_j;
    b_j   = b_jm1;
  }
  
  *b_n   = b_jp1;
  *b_nm1 = b_j;
}

/* checked OK [GJ] Sun May  3 15:43:00 EDT 1998 */
static inline int taylor_In_scaled_4(const int n, const double x, double * b_n, double * b_nm1)
{
  double ex = exp(-x);
  int status1 = gsl_sf_bessel_Inu_Jnu_taylor_impl(n, x, 1, 4, b_n);
  int status2 = GSL_SUCCESS;
  *b_n *= ex;
  if(b_nm1 != (double *)0) {
    status2 = gsl_sf_bessel_Inu_Jnu_taylor_impl(n-1, x, 1, 4, b_nm1);
    *b_nm1 *= ex;
  }
  if(status1 == GSL_EUNDRFLW || status2 == GSL_EUNDRFLW) return GSL_EUNDRFLW;
  else return GSL_SUCCESS;
}

/* convenience function for I_n and I_{n-1}
 * assumes n >= 2; x > 0
 * checked OK [GJ] Sun May  3 16:19:10 EDT 1998 
 */
static int bessel_In_scaled(const int n, const double x, double * b_n, double * b_nm1)
{
  *b_n = 0.;
  if(b_nm1 != (double *)0) *b_nm1 = 0.;

  if(x*x < 10.*(n+1)*GSL_ROOT5_MACH_EPS) {
    return taylor_In_scaled_4(n, x, b_n, b_nm1);
  }
  else if( Min( 0.29/(n*n), 0.5/(n*n + x*x) ) < GSL_ROOT3_MACH_EPS) {
    gsl_sf_bessel_Inu_scaled_asymp_unif_impl(n, x, b_n);
    if(b_nm1 != (double *)0) {
      gsl_sf_bessel_Inu_scaled_asymp_unif_impl(n-1, x, b_nm1);
    }
    return GSL_SUCCESS;
  }
  else {
    double local_b_n;
    double local_b_nm1;
    asymp_recurse_In_scaled(n, x, &local_b_n, &local_b_nm1);
    *b_n = local_b_n;
    if(b_nm1 != (double *) 0) *b_nm1 = local_b_nm1;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* checked OK [GJ] Sun May  3 15:50:23 EDT 1998 */
int gsl_sf_bessel_In_scaled_impl(int n, const double x, double * result)
{
  n = abs(n);  /* I(-n, z) = I(n, z) */
 
  if(n == 0) {
    return gsl_sf_bessel_I0_scaled_impl(x, result);
  }
  else if(n == 1) {
    return gsl_sf_bessel_I1_scaled_impl(x, result);
  }
  else { 
    if(x == 0.) {
      *result = 0.;
      return GSL_SUCCESS;
    }
    else {
      int status = bessel_In_scaled(n, fabs(x), result, (double *)0);
      if(x < 0. && GSL_IS_ODD(n)) *result = - *result;
      return status;
    }
  }
}

/* checked OK [GJ] Sun May  3 16:04:13 EDT 1998 */
int gsl_sf_bessel_In_scaled_array_impl(int n, const double x, double * result_array)
{
  n = abs(n);  /* I(-n, z) = I(n, z) */
 
  if(n == 0) {
    return gsl_sf_bessel_I0_scaled_impl(x, &(result_array[0]));
  }
  else if(n == 1) {
    int status = gsl_sf_bessel_I1_scaled_impl(x, &(result_array[1]));
    gsl_sf_bessel_I0_scaled_impl(x, &(result_array[0]));
    return status;
  }
  else { 
    if(x == 0.) {
      int j;
      result_array[0] = 1.;
      for(j=1; j<=n; j++) result_array[j] = 0.;
      return GSL_SUCCESS;
    }
    else {
      int j;
      double b_n, b_nm1;
      double ax = fabs(x);
      double two_over_x = 2./ax;
      double b_jm1, b_j, b_jp1;
      bessel_In_scaled(n, ax, &b_n, &b_nm1);  /* starting values */

      result_array[n] = b_n;
      b_jp1 = b_n;
      b_j   = b_nm1;
      for(j=n-1; j>=0; j--){
        b_jm1 = b_jp1 + j * two_over_x * b_j;
        b_jp1 = b_j;
        b_j   = b_jm1;
        result_array[j] = b_jp1;
      }
      
      /* deal with signs */
      if(x < 0.) {
        for(j=n; j>=0; j--) {
	  if(GSL_IS_ODD(j)) result_array[j] = -result_array[j];
	}
      }

      return GSL_SUCCESS;
    }
  }
}

/* checked OK [GJ] Sun May  3 16:18:04 EDT 1998 */
int gsl_sf_bessel_In_impl(int n, const double x, double * result)
{
  int j;
  double ax = fabs(x);
  
  n = abs(n);  /* I(-n, z) = I(n, z) */
  
  if(ax > GSL_LOG_DBL_MAX - 1.) {
    *result = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    int status = bessel_In_scaled(n, ax, result, (double *)0);
    if(x < 0. && GSL_IS_ODD(n)) *result = - *result;
    *result *= exp(ax);
    return status;
  }
}

/* checked OK [GJ] Sun May  3 16:08:25 EDT 1998 */
int gsl_sf_bessel_In_array_impl(int n, const double x, double * result_array)
{
  int j;
  double ax = fabs(x);
  
  n = abs(n);  /* I(-n, z) = I(n, z) */

  if(ax > GSL_LOG_DBL_MAX - 1.) {
    for(j=0; j<=n; j++) result_array[j] = 0.; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    double y;
    double eax = exp(ax);
    int status = gsl_sf_bessel_In_scaled_array_impl(n, x, result_array);
    for(j=0; j<=n; j++) result_array[j] *= eax;
    return status;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_In_scaled_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_In_scaled_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_In_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_e(const int n, const double x, double * result)
{
  int status = gsl_sf_bessel_In_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_In_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_scaled_array_e(const int n, const double x, double * result_array)
{
  int status = gsl_sf_bessel_In_scaled_array_impl(n, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_In_scaled_array_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_array_e(const int n, const double x, double * result_array)
{
  int status = gsl_sf_bessel_In_impl(n, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_In_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_In_scaled(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_In_scaled_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_In_scaled", status);
  }
  return y;
}

double gsl_sf_bessel_In(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_In_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_In", status);
  }
  return y;
}
