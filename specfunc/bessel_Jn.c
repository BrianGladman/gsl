/* Author:  G. Jungman
 * RCS:     $Id$
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The ratio of jn(x) to j0(x) is computed by backward
 * recurrence.  First the ratio jn/jn-1 is found by a
 * continued fraction expansion.  Then the recurrence
 * relating successive orders is applied until j0 or j1 is
 * reached.
 *

 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   range      # trials      peak         rms
 *    DEC       0, 30        5500       6.9e-17     9.3e-18
 *    IEEE      0, 30        5000       4.4e-16     7.9e-17
 *
 *
 * Not suitable for large n or x. Use jv() instead.
 *
 */

/*							jn.c
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_bessel.h"

extern int gsl_sf_bessel_J0_impl(const double, double *);
extern int gsl_sf_bessel_J1_impl(const double, double *);
extern int gsl_sf_fact_impl(const int, double *);


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_bessel_Jn_impl(int n, double x, double * result)
{
  int sign = 1;

  if(n < 0) {
    /* reduce to case n >= 0 */
    n = -n;
    if(GSL_IS_ODD(n)) sign = -1;
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
  else if(n == 2) {
    if(x == 0.) {
      *result = 0.;
      return GSL_SUCCESS;
    }
    else if(x < GSL_SQRT_MACH_EPS) {
      *result = 0.25*x*x;
      if(*result == 0.) {
        return GSL_EUNDRFLW;
      }
      else {
	return GSL_SUCCESS;
      }  
    }
    else {
      double b0, b1;
      gsl_sf_bessel_J0_impl(x, &b0);
      gsl_sf_bessel_J1_impl(x, &b1);
      *result = sign * (2.0 * b1 / x  -  b0);
      return GSL_SUCCESS;
    }    
  }
  else {
    if(x == 0.) {
      *result = 0.;
      return GSL_SUCCESS;
    }
    else if(x < 2.*GSL_SQRT_MACH_EPS * n) {
      double nfact;
      gsl_sf_fact_impl(n, &nfact);
      *result = gsl_sf_pow_int(0.5*x, n)/nfact;
      if(*result == 0.) {
	return GSL_EUNDRFLW;
      }
      else {
	return GSL_SUCCESS;
      }
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
	ans = gsl_sf_bessel_J1(x)/pk;
      }
      else {
	ans = gsl_sf_bessel_J0(x)/pkm1;
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
    GSL_WARNING("gsl_sf_bessel_Jn");
  }
  return y;
}
