/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "bessel.h"
#include "gsl_sf_bessel.h"


/* I_{n+1}/I_n
 */
static
int
bessel_In_CF1(const int n, const double x, const double threshold, double * ratio)
{
  const int kmax = 2000;
  double tk   = 1.0;
  double sum  = 1.0;
  double rhok = 0.0;
  int k;

  for(k=1; k<=kmax; k++) {
    double ak = 0.25 * (x/(n+k)) * (x/(n+k+1.0));
    rhok = -ak*(1.0 + rhok)/(1.0 + ak*(1.0 + rhok));
    tk  *= rhok;
    sum += tk;
    if(fabs(tk/sum) < threshold) break;
  }

  *ratio = 0.5*x/(n+1.0) * sum;

  if(k == kmax)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_bessel_In_scaled_impl(int n, const double x, gsl_sf_result * result)
{
  const double ax = fabs(x);

  n = abs(n);  /* I(-n, z) = I(n, z) */

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(n == 0) {
    return gsl_sf_bessel_I0_scaled_impl(x, result);
  }
  else if(n == 1) {
    return gsl_sf_bessel_I1_scaled_impl(x, result);
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x*x < 10.0*(n+1.0)/M_E) {
    double t;
    int stat_In = gsl_sf_bessel_Inu_Jnu_taylor_impl((double)n, ax, 1, 50, GSL_DBL_EPSILON, &t);
    result->val = t * exp(-ax);
    result->err = GSL_DBL_EPSILON * result->val * fabs(ax);
    if(x < 0.0 && GSL_IS_ODD(n)) result->val = -result->val;
    return stat_In;
  }
  else if(n < 150) {
    gsl_sf_result I0_scaled;
    int stat_I0 = gsl_sf_bessel_I0_scaled_impl(ax, &I0_scaled);
    double rat;
    int stat_CF1 = bessel_In_CF1(n, ax, GSL_DBL_EPSILON, &rat);
    double Ikp1 = rat * GSL_SQRT_DBL_MIN;
    double Ik	= GSL_SQRT_DBL_MIN;
    double Ikm1;
    int k;
    for(k=n; k >= 1; k--) {
      Ikm1 = Ikp1 + 2.0*k/ax * Ik;
      Ikp1 = Ik;
      Ik   = Ikm1;
    }
    result->val = I0_scaled.val * (GSL_SQRT_DBL_MIN / Ik);
    result->err = I0_scaled.err * (GSL_SQRT_DBL_MIN / Ik);
    if(x < 0.0 && GSL_IS_ODD(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_2(stat_I0, stat_CF1);
  }
  else if( GSL_MIN( 0.29/(n*n), 0.5/(n*n + x*x) ) < 0.5*GSL_ROOT3_DBL_EPSILON) {
    int stat_as = gsl_sf_bessel_Inu_scaled_asymp_unif_impl((double)n, ax, result);
    if(x < 0.0 && GSL_IS_ODD(n)) result->val = -result->val;
    return stat_as;
  }
  else {
    const int nhi = 2 + (int) (1.2 / GSL_ROOT6_DBL_EPSILON);
    gsl_sf_result r_Ikp1;
    gsl_sf_result r_Ik;
    int stat_a1 = gsl_sf_bessel_Inu_scaled_asymp_unif_impl(nhi+1.0,     ax, &r_Ikp1);
    int stat_a2 = gsl_sf_bessel_Inu_scaled_asymp_unif_impl((double)nhi, ax, &r_Ik);
    double Ikp1 = r_Ikp1.val;
    double Ik   = r_Ik.val;
    double Ikm1;
    int k;
    for(k=nhi; k > n; k--) {
      Ikm1 = Ikp1 + 2.0*k/ax * Ik;
      Ikp1 = Ik;
      Ik   = Ikm1;
    }
    result->val = Ik;
    result->err = Ik * (r_Ikp1.err/r_Ikp1.val + r_Ik.err/r_Ik.val);
    if(x < 0.0 && GSL_IS_ODD(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_2(stat_a1, stat_a2);
  }
}


int
gsl_sf_bessel_In_scaled_array_impl(const int nmin, const int nmax, const double x, double * result_array)
{
  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(nmax < nmin || nmin < 0) {
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
    gsl_sf_result I0_scaled;
    int stat = gsl_sf_bessel_I0_scaled_impl(x, &I0_scaled);
    result_array[0] = I0_scaled.val;
    return stat;
  }
  else {
    const double ax = fabs(x);
    const double two_over_x = 2.0/ax;

    /* starting values */
    gsl_sf_result r_Inp1;
    gsl_sf_result r_In;
    int stat_0 = gsl_sf_bessel_In_scaled_impl(nmax+1, ax, &r_Inp1);
    int stat_1 = gsl_sf_bessel_In_scaled_impl(nmax,   ax, &r_In);
    double Inp1 = r_Inp1.val;
    double In   = r_In.val;
    double Inm1;
    int n;

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

    return GSL_ERROR_SELECT_2(stat_0, stat_1);
  }
}


int
gsl_sf_bessel_In_impl(const int n_in, const double x, gsl_sf_result * result)
{
  const double ax = fabs(x);
  const int n = abs(n_in);  /* I(-n, z) = I(n, z) */
  gsl_sf_result In_scaled;
  const int stat_In_scaled = gsl_sf_bessel_In_scaled_impl(n, ax, &In_scaled);

  /* In_scaled is always less than 1,
   * so this overflow check is conservative.
   */
  if(ax > GSL_LOG_DBL_MAX - 1.0) {
    result->val = 0.0; /* FIXME: should be Inf */
    result->err = 0.0;
    return GSL_EOVRFLW;
  }
  else {
    const double ex = exp(ax);
    result->val = ex * In_scaled.val;
    result->err = ex * In_scaled.err;
    if(x < 0.0 && GSL_IS_ODD(n)) result->val = -result->val;
    return stat_In_scaled;
  }
}


int
gsl_sf_bessel_In_array_impl(const int nmin, const int nmax, const double x, double * result_array)
{
  double ax = fabs(x);

  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(ax > GSL_LOG_DBL_MAX - 1.0) {
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

int gsl_sf_bessel_In_scaled_e(const int n, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_In_scaled_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_In_scaled_e", status);
  }
  return status;
}

int gsl_sf_bessel_In_e(const int n, const double x, gsl_sf_result * result)
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
