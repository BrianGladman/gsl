/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_bessel.h"


double gsl_sf_bessel_I_scaled(int n, double x)
{
  double result;
  n = abs(n);  /* I(-n, z) = I(n, z) */
 
  if(n == 0) {
    result = gsl_sf_bessel_I0_scaled(x);
  }
  else if(n == 1) {
    result = gsl_sf_bessel_I1_scaled(x);
  }
  else { 
    if(x == 0.) {
      result = 0.;
    }
    else {
      int j;
      int jmax = 2 * (n + (int) sqrt(50.*n));
      double renorm     = 1.e+12;
      double renorm_inv = 1.e-12;
      double two_over_x = 2./fabs(x);
      double b_np1 = 0.;
      double b_jp1 = 0.;
      double b_j   = 1.;
      double b_jm1;
    
      /* Downward recursion. [Gradshteyn + Ryzhik, 8.471.1] */
      for(j=jmax; j>0; j--){
      	b_jm1 = b_jp1 + j * two_over_x * b_j;
      	b_jp1 = b_j;
      	b_j   = b_jm1;
      
      	/* Renormalize to prevent overflow. */
      	if(fabs(b_j) > renorm){
	  b_np1 *= renorm_inv;
	  b_j   *= renorm_inv;
	  b_jp1 *= renorm_inv;
      	}
      
	/* Grab the value that we want on the way down. */
      	if(j == n) b_np1 = b_jp1;
      }

      /* Normalize using I0(x); j=0 here. */
      b_np1 *= gsl_sf_bessel_I0_scaled(x) / b_j;

      result = x < 0. && (GSL_IS_ODD(n)) ? -b_np1 : b_np1;
    }
  }
  
  return result;
}


double gsl_sf_bessel_I(int n, double x)
{
  double y = fabs(x);
  
  if(y > GSL_LOG_DBL_MAX) {
    /* domain error */
    return 0.;
  }
  else {
    return exp(y) * gsl_sf_bessel_I_scaled(n, x);
  }
}
