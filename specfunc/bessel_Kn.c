/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_besel_.h"


double gsl_sf_bessel_K_scaled(int n, double x)
{
  double result;
  
  if(x <= 0.) {
    /* domain error */
    result = 0.;
  }
  else {
    n = abs(n); /* K(-n, z) = K(n, z) */
    
    if(n == 0) {
      result = gsl_sf_bessel_K0_scaled(x);
    }
    else if(n == 1) {
      result = gsl_sf_bessel_K1_scaled(x);
    }
    else {
      /* Upward recurrence. [Gradshteyn + Ryzhik, 8.471.1] */
      int j;
      double two_over_x = 2.0/x;
      double b_jm1 = gsl_sf_bessel_K0_scaled(x);
      double b_j   = gsl_sf_bessel_K1_scaled(x);
      double b_jp1;

      for(j=1; j<n; j++) {
	b_jp1 = b_jm1 + j * two_over_x * b_j;
	b_jm1 = b_j;
	b_j   = b_jp1; 
      } 
      
      result = b_j; 
    }
  }
  
  return result;
}


double gsl_sf_bessel_K(int n, double x)
{
  if(x <= 0.) {
    /* domain error */
    return 0.;
  }
  else {
    return exp(-x) * gsl_sf_bessel_K_scaled(n, x);
  }
}
