/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* backward recursion using asymptotic starting values
 * assumes: n >= 2, x > 0
 * checked OK [GJ] Sun May  3 09:15:06 EDT 1998 
 */
static
void
asymp_recurse_In_scaled(const int n, const double x, double * b_n, double * b_nm1)
{
  int j;
  double b_jp1, b_j, b_jm1;
  double two_over_x = 2./x;
  
  /* go high enough to apply asymptotic result */
  int big_nu = 2 + GSL_MAX_INT(n, (int)sqrt(GSL_MAX_DBL(0.0, 0.5/GSL_ROOT3_MACH_EPS - x*x)));

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
#ifdef HAVE_INLINE
inline
#endif
static
int
taylor_In_scaled(const int max, const int n, const double x, double * b_n, double * b_nm1)
{
  double ex = exp(-x);
  int status1 = gsl_sf_bessel_Inu_Jnu_taylor_impl(n, x, 1, max, b_n);
  int status2 = GSL_SUCCESS;
  *b_n *= ex;
  if(b_nm1 != (double *)0) {
    status2 = gsl_sf_bessel_Inu_Jnu_taylor_impl(n-1, x, 1, max, b_nm1);
    *b_nm1 *= ex;
  }
  return GSL_ERROR_SELECT_2(status1, status2);
}


/* convenience function for I_n and I_{n-1}
 * assumes n >= 2; x > 0
 * checked OK [GJ] Sun May  3 16:19:10 EDT 1998 
 */
static
int
bessel_In_scaled(const int n, const double x, double * b_n, double * b_nm1)
{
  *b_n = 0.0;
  if(b_nm1 != (double *)0) *b_nm1 = 0.0;

  if(x*x < 10.0*(n+1)*GSL_ROOT5_MACH_EPS) {
    return taylor_In_scaled(4, n, x, b_n, b_nm1);
  }
  else if(x*x < 10.0*(n+1)/M_E) {
    return taylor_In_scaled(14, n, x, b_n, b_nm1);
  }
  else if( GSL_MIN( 0.29/(n*n), 0.5/(n*n + x*x) ) < GSL_ROOT3_MACH_EPS) {
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
int
gsl_sf_bessel_In_scaled_impl(int n, const double x, double * result)
{
  n = abs(n);  /* I(-n, z) = I(n, z) */
 
  if(n == 0) {
    return gsl_sf_bessel_I0_scaled_impl(x, result);
  }
  else if(n == 1) {
    return gsl_sf_bessel_I1_scaled_impl(x, result);
  }
  else { 
    if(x == 0.0) {
      *result = 0.0;
      return GSL_SUCCESS;
    }
    else {
      int status = bessel_In_scaled(n, fabs(x), result, (double *)0);
      if(x < 0.0 && GSL_IS_ODD(n)) *result = - *result;
      return status;
    }
  }
}


/* checked OK [GJ] Sun May  3 16:04:13 EDT 1998 */
int
gsl_sf_bessel_In_scaled_array_impl(const int nmin, const int nmax, const double x, double * result_array)
{
  if(nmax < nmin || nmin < 0) {
    int j;
    for(j=0; j<=nmax-nmin; j++) result_array[j] = 0.0;
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    int j;
    for(j=0; j<=nmax-nmin; j++) result_array[j] = 0.0;
    if(nmin == 0) result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(nmax == 0) {
    return gsl_sf_bessel_I0_scaled_impl(x, &(result_array[0]));
  }
  else {
    const double ax = fabs(x);
    const double two_over_x = 2.0/ax;
    double Inm1;
    double In;
    double Inp1;
    int n;

    int stat = bessel_In_scaled(nmax+1, ax, &Inp1, &In);  /* starting values */

    for(n=nmax; n>=nmin; n--) {
      result_array[n-nmin] = In;
      Inm1 = Inp1 + n * two_over_x * In;
      Inp1 = In;
      In   = Inm1;
    }

    /* deal with signs */
    if(x < 0.0) {
      for(n=nmin; n<=nmax; n++) {
        if(GSL_IS_ODD(n)) result_array[n-nmin] = -result_array[n-nmin];
      }
    }

    return stat;
  }
}

/* checked OK [GJ] Sun May  3 16:18:04 EDT 1998 */
int
gsl_sf_bessel_In_impl(int n, const double x, double * result)
{
  double ax = fabs(x);
  
  n = abs(n);  /* I(-n, z) = I(n, z) */
  
  if(ax > GSL_LOG_DBL_MAX - 1.0) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    int status = bessel_In_scaled(n, ax, result, (double *)0);
    if(x < 0.0 && GSL_IS_ODD(n)) *result = - *result;
    *result *= exp(ax);
    return status;
  }
}

/* checked OK [GJ] Sun May  3 16:08:25 EDT 1998 */
int
gsl_sf_bessel_In_array_impl(const int nmin, const int nmax, const double x, double * result_array)
{
  double ax = fabs(x);

  if(ax > GSL_LOG_DBL_MAX - 1.0) {
    int j;
    for(j=0; j<=nmax-nmin; j++) result_array[j] = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else {
    int j;
    double eax = exp(ax);
    int status = gsl_sf_bessel_In_scaled_array_impl(nmin, nmax, x, result_array);
    for(j=0; j<=nmax-nmin; j++) result_array[j] *= eax;
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

int
gsl_sf_bessel_In_scaled_array_e(const int nmin, const int nmax, const double x, double * result_array)
{
  int status = gsl_sf_bessel_In_scaled_array_impl(nmin, nmax, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_In_scaled_array_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_array_e(const int nmin, const int nmax, const double x, double * result_array)
{
  int status = gsl_sf_bessel_In_array_impl(nmin, nmax, x, result_array);
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
