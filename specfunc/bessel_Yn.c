/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "bessel_amp_phase.h"
#include "bessel_olver.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_psi.h"
#include "gsl_sf_bessel.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* assumes n >= 1 */
static int bessel_Yn_small_x(const int n, const double x, gsl_sf_result * result)
{
  int k;
  double y = 0.25 * x * x;
  double ln_x_2 = log(0.5*x);
  gsl_sf_result ln_nm1_fact;
  double k_term;
  double term1, sum1, ln_pre1;
  double term2, sum2, pre2;

  gsl_sf_lnfact_impl((unsigned int)(n-1), &ln_nm1_fact);

  ln_pre1 = -n*ln_x_2 + ln_nm1_fact.val;
  if(ln_pre1 > GSL_LOG_DBL_MAX - 3.0) return GSL_EOVRFLW;

  sum1 = 1.0;
  k_term = 1.0;
  for(k=1; k<=n-1; k++) {
    k_term *= y/(k * (n-k));
    sum1 += k_term;
  }
  term1 = -exp(ln_pre1) * sum1 / M_PI;
  
  pre2 = -exp(n*ln_x_2) / M_PI;
  if(fabs(pre2) > 0.0) {
    const int KMAX = 20;
    gsl_sf_result psi_n;
    gsl_sf_result npk_fact;
    double yk = 1.0;
    double k_fact  = 1.0;
    double psi_kp1 = -M_EULER;
    double psi_npkp1;
    gsl_sf_psi_int_impl(n, &psi_n);
    gsl_sf_fact_impl((unsigned int)n, &npk_fact);
    psi_npkp1 = psi_n.val + 1.0/n;
    sum2 = (psi_kp1 + psi_npkp1 - 2.0*ln_x_2)/npk_fact.val;
    for(k=1; k<KMAX; k++) {
      psi_kp1   += 1./k;
      psi_npkp1 += 1./(n+k);
      k_fact   *= k;
      npk_fact.val *= n+k;
      yk *= -y;
      k_term = yk*(psi_kp1 + psi_npkp1 - 2.0*ln_x_2)/(k_fact*npk_fact.val);
      sum2 += k_term;
    }
    term2 = pre2 * sum2;
  }
  else {
    term2 = 0.0;
  }

  result->val  = term1 + term2;
  result->err  = GSL_DBL_EPSILON * (fabs(ln_pre1)*fabs(term1) + fabs(term2));
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_bessel_Yn_impl(int n, const double x, gsl_sf_result * result)
{
  int sign = 1;

  if(n < 0) {
    /* reduce to case n >= 0 */
    n = -n;
    if(GSL_IS_ODD(n)) sign = -1;
  }

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(n == 0) {
    int status = gsl_sf_bessel_Y0_impl(x, result);
    result->val *= sign;
    return status;
  }
  else if(n == 1) {
    int status = gsl_sf_bessel_Y1_impl(x, result);
    result->val *= sign;
    return status;
  }
  else {
    if(x <= 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_EDOM;
    }
    if(x < 5.0) {
      int status = bessel_Yn_small_x(n, x, result);
      result->val *= sign;
      return status;
    }
    else if(GSL_ROOT3_DBL_EPSILON * x > (n*n + 1.0)) {
      int status = gsl_sf_bessel_Ynu_asympx_impl((double)n, x, result);
      result->val *= sign;
      return status;
    }
    else if(n > 50) {
      int status = gsl_sf_bessel_Ynu_asymp_Olver_impl((double)n, x, result);
      result->val *= sign;
      return status;
    }
    else {
      double two_over_x = 2.0/x;
      gsl_sf_result r_by;
      gsl_sf_result r_bym;
      int stat_1 = gsl_sf_bessel_Y1_impl(x, &r_by);
      int stat_0 = gsl_sf_bessel_Y0_impl(x, &r_bym);
      double bym = r_bym.val;
      double by  = r_by.val;
      double byp;
      int j;

      for(j=1; j<n; j++) { 
	byp = j*two_over_x*by - bym;
	bym = by;
	by  = byp;
      }
      result->val  = sign * by;
      result->err  = fabs(result->val) * (fabs(r_by.err/r_by.val) + fabs(r_bym.err/r_bym.val));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_ERROR_SELECT_2(stat_1, stat_0);
    }
  }
}


int
gsl_sf_bessel_Yn_array_impl(const int nmin, const int nmax, const double x, double * result_array)
{
  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(nmin < 0 || nmax < nmin || x <= 0.0) {
    int j;
    for(j=0; j<=nmax-nmin; j++) result_array[j] = 0.0;
    return GSL_EDOM;
  }
  else {
    gsl_sf_result r_Ynm1;
    gsl_sf_result r_Yn;
    int stat_nm1 = gsl_sf_bessel_Yn_impl(nmin,   x, &r_Ynm1);
    int stat_n   = gsl_sf_bessel_Yn_impl(nmin+1, x, &r_Yn);
    double Ynp1;
    double Yn   = r_Yn.val;
    double Ynm1 = r_Ynm1.val;
    int n;

    int stat = GSL_ERROR_SELECT_2(stat_nm1, stat_n);

    if(stat == GSL_SUCCESS) {
      for(n=nmin+1; n<=nmax+1; n++) {
        result_array[n-nmin-1] = Ynm1;
        Ynp1 = -Ynm1 + 2.0*n/x * Yn;
	Ynm1 = Yn;
        Yn   = Ynp1;
      }
    }
    else {
      for(n=nmin; n<=nmax; n++) {
        result_array[n-nmin] = 0.0;
      }
    }

    return stat;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_bessel_Yn_e(const int n, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Yn_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Yn_e", status);
  }
  return status;
}


int
gsl_sf_bessel_Yn_array_e(const int nmin, const int nmax, const double x, double * result_array)
{
  int status = gsl_sf_bessel_Yn_array_impl(nmin, nmax, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Yn_array_e", status);
  }
  return status;
}
