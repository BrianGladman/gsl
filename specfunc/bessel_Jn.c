/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "bessel.h"
#include "bessel_amp_phase.h"
#include "bessel_olver.h"
#include "gsl_sf_bessel.h"



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int gsl_sf_bessel_Jn_impl(int n, double x, gsl_sf_result * result)
{
  int sign = 1;

  if(n < 0) {
    /* reduce to case n >= 0 */
    n = -n;
    if(GSL_IS_ODD(n)) sign = -sign;
  }  

  if(x < 0.0) {
    /* reduce to case x >= 0. */
    x = -x;
    if(GSL_IS_ODD(n)) sign = -sign;
  }

  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(n == 0) {
    gsl_sf_result b0;
    int stat_J0 = gsl_sf_bessel_J0_impl(x, &b0);
    result->val = sign * b0.val;
    result->err = b0.err;
    return stat_J0;
  }
  else if(n == 1) {
    gsl_sf_result b1;
    int stat_J1 = gsl_sf_bessel_J1_impl(x, &b1);
    result->val = sign * b1.val;
    result->err = b1.err;
    return stat_J1;
  }
  else {
    if(x == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(x*x < 10.0*(n+1.0)*GSL_ROOT5_DBL_EPSILON) {
      gsl_sf_result b;
      int status = gsl_sf_bessel_IJ_taylor_impl((double)n, x, -1, 50, GSL_DBL_EPSILON, &b);
      result->val  = sign * b.val;
      result->err  = b.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return status;
    }
    else if(GSL_ROOT3_DBL_EPSILON * x > (n*n+1.0)) {
      int status = gsl_sf_bessel_Jnu_asympx_impl((double)n, x, result);
      result->val *= sign;
      return status;
    }
    else if(n > 50) {
      int status = gsl_sf_bessel_Jnu_asymp_Olver_impl((double)n, x, result);
      result->val *= sign;
      return status;
    }
    else {
      double ans;
      double err;
      double ratio;
      double sgn;
      int stat_b;
      int stat_CF1 = gsl_sf_bessel_J_CF1((double)n, x, &ratio, &sgn);

      /* backward recurrence */
      double Jkp1 = GSL_SQRT_DBL_MIN * ratio;
      double Jk   = GSL_SQRT_DBL_MIN;
      double Jkm1;
      int k;

      for(k=n; k>0; k--) {
	Jkm1 = 2.0*k/x * Jk - Jkp1;
	Jkp1 = Jk;
	Jk   = Jkm1;
      }

      if(fabs(Jkp1) > fabs(Jk)) {
        gsl_sf_result b1;
	stat_b = gsl_sf_bessel_J1_impl(x, &b1);
	ans = b1.val/Jkp1 * GSL_SQRT_DBL_MIN;
	err = b1.err/Jkp1 * GSL_SQRT_DBL_MIN;
      }
      else {
        gsl_sf_result b0;
	stat_b = gsl_sf_bessel_J0_impl(x, &b0);
	ans = b0.val/Jk * GSL_SQRT_DBL_MIN;
	err = b0.err/Jk * GSL_SQRT_DBL_MIN;
      }

      result->val = sign * ans;
      result->err = fabs(err);
      return GSL_ERROR_SELECT_2(stat_CF1, stat_b);
    }
  }
}


int
gsl_sf_bessel_Jn_array_impl(int nmin, int nmax, double x, double * result_array)
{
  if(result_array == 0) {
    return GSL_EFAULT;
  }
  else if(nmin < 0 || nmax < nmin) {
    int n;
    for(n=nmax; n>=nmin; n--) {
      result_array[n-nmin] = 0.0;
    }
    return GSL_EDOM;
  }
  else if(x == 0.0) {
    int n;
    for(n=nmax; n>=nmin; n--) {
      result_array[n-nmin] = 0.0;
    }
    if(nmin == 0) result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result r_Jnp1;
    gsl_sf_result r_Jn;
    int stat_np1 = gsl_sf_bessel_Jn_impl(nmax+1, x, &r_Jnp1);
    int stat_n   = gsl_sf_bessel_Jn_impl(nmax,   x, &r_Jn);
    int stat = GSL_ERROR_SELECT_2(stat_np1, stat_n);

    double Jnp1 = r_Jnp1.val;
    double Jn   = r_Jn.val;
    double Jnm1;
    int n;

    if(stat == GSL_SUCCESS) {
      for(n=nmax; n>=nmin; n--) {
        result_array[n-nmin] = Jn;
        Jnm1 = -Jnp1 + 2.0*n/x * Jn;
        Jnp1 = Jn;
        Jn   = Jnm1;
      }
    }
    else {
      for(n=nmax; n>=nmin; n--) {
        result_array[n-nmin] = 0.0;
      }
    }

    return stat;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Jn_e(const int n, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_bessel_Jn_impl(n, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Jn_e", status);
  }
  return status;
}

int
gsl_sf_bessel_Jn_array_e(int nmin, int nmax, double x, double * result_array)
{
  int status = gsl_sf_bessel_Jn_array_impl(nmin, nmax, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_bessel_Jn_array_e", status);
  }
  return status;
}
