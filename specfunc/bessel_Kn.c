/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"

#include "psi_impl.h"
#include "gamma_impl.h"
#include "bessel.h"
#include "bessel_K0_impl.h"
#include "bessel_K1_impl.h"

#include "bessel_Kn_impl.h"

#define Min(a,b) ((a) < (b) ? (a) : (b))


/* [Abramowitz+Stegun, 9.6.11] 
 * checked OK [GJ] Sun May  3 20:48:26 EDT 1998 
 */
static int bessel_Kn_scaled_small_x(const int n, const double x, double * result)
{
  int k;
  double y = 0.25 * x * x;
  double ln_x_2 = log(0.5*x);
  double ln_nm1_fact;
  double k_term;
  double term1, sum1, ln_pre1;
  double term2, sum2, pre2;

  gsl_sf_lnfact_impl(n-1, &ln_nm1_fact);

  ln_pre1 = -n*ln_x_2 + ln_nm1_fact;
  if(ln_pre1 > GSL_LOG_DBL_MAX - 3.) return GSL_EOVRFLW;

  sum1 = 1.0;
  k_term = 1.0;
  for(k=1; k<=n-1; k++) {
    k_term *= -y/(k * (n-k));
    sum1 += k_term;
  }
  term1 = 0.5 * exp(ln_pre1) * sum1;
  
  pre2 = 0.5 * exp(n*ln_x_2);
  if(pre2 > 0.0) {
    const int KMAX = 20;
    double psi_n;
    double npk_fact;
    double yk = 1.0;
    double k_fact  = 1.0;
    double psi_kp1 = -M_EULER;
    double psi_npkp1;
    gsl_sf_psi_int_impl(n, &psi_n);
    gsl_sf_fact_impl(n, &npk_fact);
    psi_npkp1 = psi_n + 1./n;
    sum2 = (psi_kp1 + psi_npkp1 - 2.0*ln_x_2)/npk_fact;
    for(k=1; k<KMAX; k++) {
      psi_kp1   += 1./k;
      psi_npkp1 += 1./(n+k);
      k_fact   *= k;
      npk_fact *= n+k;
      yk *= y;
      k_term = yk*(psi_kp1 + psi_npkp1 - 2.0*ln_x_2)/(k_fact*npk_fact);
      sum2 += k_term;
    }
    term2 = ( GSL_IS_ODD(n) ? -1. : 1. ) * pre2 * sum2;
  }
  else {
    term2 = 0.0;
  }

  *result = exp(x) * (term1 + term2);
  return GSL_SUCCESS;
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
  else if(x <= 5.0) {
    return bessel_Kn_scaled_small_x(n, x, result);
  }
  else if(GSL_ROOT3_MACH_EPS * x > 0.25 * (n*n + 1)) {
  printf("asymp: ");
    return gsl_sf_bessel_Knu_scaled_asympx_impl((double)n, x, result);
  }
  else if(Min( 0.29/(n*n), 0.5/(n*n + x*x) ) < GSL_ROOT3_MACH_EPS) {
  printf("aunif: ");
    return gsl_sf_bessel_Knu_scaled_asymp_unif_impl((double)n, x, result);
  }
  else {
  int vvv = printf("recur: ");
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

/* checked OK [GJ] Sun May  3 21:50:17 EDT 1998 */
int gsl_sf_bessel_Kn_scaled_array_impl(int n, const double x, double * result_array)
{
  n = abs(n); /* K(-n, z) = K(n, z) */
  
  if(x <= 0.) return GSL_EDOM;

  if(n == 0) {
    return gsl_sf_bessel_K0_scaled_impl(x, &(result_array[0]));
  }
  else if(n == 1) {
    gsl_sf_bessel_K0_scaled_impl(x, &(result_array[0]));
    return gsl_sf_bessel_K1_scaled_impl(x, &(result_array[1]));
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

      result_array[j]	= b_j;
      result_array[j+1] = b_jp1;

      b_jm1 = b_j;
      b_j   = b_jp1; 
    }

    /* FIXME: this may not be portable... 
       want to test if something overflowed; note that it's enough to test
       the last element since the sequence increases monotonically
       */
    return (result_array[n] < DBL_MAX ? GSL_SUCCESS : GSL_EOVRFLW);
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
