/* Author:  G. Jungman
 * RCS:     $Id$
 *
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "bessel.h"
#include "bessel_amp_phase.h"
#include "gsl_sf_bessel.h"


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
    gsl_sf_bessel_J0_impl(x, &b0);
    *result = sign * b0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    double b1;
    gsl_sf_bessel_J1_impl(x, &b1);
    *result = sign * b1;
    return GSL_SUCCESS;
  }
  else {
    if(x == 0.) {
      *result = 0.;
      return GSL_SUCCESS;
    }
    else if(x*x < 10.*(n+1)*GSL_ROOT5_MACH_EPS) {
      int status = gsl_sf_bessel_Inu_Jnu_taylor_impl((double)n, x, -1, 4, result);
      *result *= sign;
      return status;
    }
    else if(GSL_ROOT3_MACH_EPS * x > (n*n+1)) {
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
      double pkm2, pkm1, r;

      /* continued fraction [Abramowitz+Stegun, 9.1.73] */
      int k = 60;
      double pk = 2 * (n + k);
      double xk = x * x;
      double ans = pk;

      do {
	pk -= 2.0;
	ans = pk - (xk/ans);
      }
      while(--k > 0);

      ans = x/ans;

      /* backward recurrence */
      pk = 1.0;
      pkm1 = 1.0/ans;
      k = n-1;
      r = 2 * k;

      do {
	pkm2 = (pkm1 * r  -  pk * x) / x;
	pk = pkm1;
	pkm1 = pkm2;
	r -= 2.0;
      }
      while( --k > 0 );

      if(fabs(pk) > fabs(pkm1)) {
        double b1;
	gsl_sf_bessel_J1_impl(x, &b1);
	ans = b1/pk;
      }
      else {
        double b0;
	gsl_sf_bessel_J0_impl(x, &b0);
	ans = b0/pkm1;
      }

      *result = sign * ans;
      return GSL_SUCCESS;
    }
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
