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


double gsl_sf_bessel_Jn(int n, double x )
{
  double result;
  int sign = 1;

  if(n < 0) {
    /* reduce to case n >= 0 */
    n = -n;
    if(GSL_IS_ODD(n)) sign = -1;
  }  

  if(x < 0.) {
    /* reduce to case x >= 0. */
    if(GSL_IS_ODD(n)) sign = -sign;
    x = -x;
  }

  if(n == 0) {
    result = sign * gsl_sf_bessel_J0(x);
  }
  else if(n == 1) {
    result = sign * gsl_sf_bessel_J1(x) );
  }
  else if(n == 2) {
    if(x < GSL_SQRT_MACH_EPS) {
      /* underflow */
      result = 0.;
    }
    else {
      result = sign * (2.0 * gsl_sf_bessel_J1(x) / x  -  gsl_sf_bessel_J0(x));
    }    
  }
  else {  
    if(x < GSl_MACH_EPS) {
      /* underflow */
      result = 0.;
    }
    else {
      double pkm2, pkm1, r;

      /* continued fraction */
      int k = 53;
      double pk = 2 * (n + k);
      double xk = x * x;
      double result = pk;

      do {
	pk -= 2.0;
	result = pk - (xk/result);
      }
      while( --k > 0 );

      result = x/result;

      /* backward recurrence */
      pk = 1.0;
      pkm1 = 1.0/result;
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
	result = gsl_sf_bessel_J1(x)/pk;
      }
      else {
	result = gsl_sf_bessel_J0(x)/pkm1;
      }

      result = sign * result;
    }
  }
}
