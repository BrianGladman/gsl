/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "bessel.h"
#include "bessel_amp_phase.h"
#include "bessel_olver.h"
#include "gsl_sf_bessel.h"

/* continued fraction [Abramowitz+Stegun, 9.1.73]
 */
static
int
bessel_Jn_CF1(const int n, const double x, double * ratio)
{
  int k = 60;
  double pk  = 2 * (n + k);
  double xk  = x * x;
  double ans = pk;

  do {
    pk -= 2.0;
    if(ans == 0.0) {
      *ratio = 0.0;
      return GSL_EZERODIV;
    }
    ans = pk - (xk/ans);
  }
  while(--k > 0);

  if(ans != 0.0) {
    *ratio = x/ans;
    return GSL_SUCCESS;
  }
  else {
    *ratio = 0.0;
    return GSL_EZERODIV;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

/* checked OK [GJ] Sun May  3 22:35:34 EDT 1998 */
int gsl_sf_bessel_Jn_impl(int n, double x, double * result)
{
  int sign = 1;

  if(n < 0) {
    /* reduce to case n >= 0 */
    n = -n;
    if(GSL_IS_ODD(n)) sign = -sign;
  }  

  if(x < 0.) {
    /* reduce to case x >= 0. */
    x = -x;
    if(GSL_IS_ODD(n)) sign = -sign;
  }

  if(n == 0) {
    double b0;
    int stat_J0 = gsl_sf_bessel_J0_impl(x, &b0);
    *result = sign * b0;
    return stat_J0;
  }
  else if(n == 1) {
    double b1;
    int stat_J1 = gsl_sf_bessel_J1_impl(x, &b1);
    *result = sign * b1;
    return stat_J1;
  }
  else {
    if(x == 0.0) {
      *result = 0.0;
      return GSL_SUCCESS;
    }
    else if(x*x < 10.0*(n+1.0)*GSL_ROOT5_MACH_EPS) {
      int status = gsl_sf_bessel_Inu_Jnu_taylor_impl((double)n, x, -1, 4, result);
      *result *= sign;
      return status;
    }
    else if(GSL_ROOT3_MACH_EPS * x > (n*n+1.0)) {
      int status = gsl_sf_bessel_Jnu_asympx_impl((double)n, x, result);
      *result *= sign;
      return status;
    }
    else if(n > 30) {
      int status = gsl_sf_bessel_Jnu_asymp_Olver_impl((double)n, x, result);
      *result *= sign;
      return status;
    }
    else {
      double ans;
      double ratio;
      int stat_b;
      int stat_CF1 = bessel_Jn_CF1(n, x, &ratio);

      /* backward recurrence */
      double Jk   = ratio;
      double Jkm1 = 1.0;
      double Jkm2;
      int k;

      for(k=n-1; k>0; k--) {
        double r = 2.0 * k;
	Jkm2 = (Jkm1 * r  -  Jk * x) / x;
	Jk   = Jkm1;
	Jkm1 = Jkm2;
      }

      if(fabs(Jk) > fabs(Jkm1)) {
        double b1;
	stat_b = gsl_sf_bessel_J1_impl(x, &b1);
	ans = b1/Jk * ratio;
      }
      else {
        double b0;
	stat_b = gsl_sf_bessel_J0_impl(x, &b0);
	ans = b0/Jkm1 * ratio;
      }

      *result = sign * ans;
      return GSL_ERROR_SELECT_2(stat_CF1, stat_b);
    }
  }
}


int
gsl_sf_bessel_Jn_array_impl(int nmin, int nmax, double x, double * result_array)
{
  if(nmin < 0 || nmax < nmin) {
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
    double Jnp1;
    double Jn;
    double Jnm1;
    int n;

    int stat_np1 = gsl_sf_bessel_Jn_impl(nmax+1, x, &Jnp1);
    int stat_n   = gsl_sf_bessel_Jn_impl(nmax,   x, &Jn);

    int stat = GSL_ERROR_SELECT_2(stat_np1, stat_n);

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

int gsl_sf_bessel_Jn_e(const int n, const double x, double * result)
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


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_bessel_Jn(const int n, const double x)
{
  double y;
  int status = gsl_sf_bessel_Jn_impl(n, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_bessel_Jn", status);
  }
  return y;
}
