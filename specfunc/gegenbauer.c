/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gegenbauer.h"


/* See: [Thompson, Atlas for Computing Mathematical Functions] */


double
gsl_sf_gegenpoly_1(double lambda, double x)
{
  if(lambda == 0.0)
    return 2.0*x;
  else
    return 2.0*lambda*x;
}

double
gsl_sf_gegenpoly_2(double lambda, double x)
{
  if(lambda == 0.0)
    return -1.0 + 2.0*x*x;
  else
    return lambda*(-1.0 + 2.0*(1.0+lambda)*x*x);
}

double
gsl_sf_gegenpoly_3(double lambda, double x)
{
  if(lambda == 0.0)
    return x*(-2.0 + 4.0/3.0*x*x);
  else {
    double c = 4.0 + lambda*(6.0 + 2.0*lambda);
    return 2.0*lambda * x * ( -1.0 - lambda + c*x*x/3.0 );
  }
}


int
gsl_sf_gegenpoly_n_impl(int n, double lambda, double x, double * result)
{
  if(lambda <= -0.5 || n < 0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  
  if(n == 0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    *result = gsl_sf_gegenpoly_1(lambda, x);
    return GSL_SUCCESS;
  }
  else if(n == 2) {
    *result = gsl_sf_gegenpoly_2(lambda, x);
    return GSL_SUCCESS;
  }
  else if(n == 3) {
    *result = gsl_sf_gegenpoly_3(lambda, x);
    return GSL_SUCCESS;
  }
  else {
    if(lambda == 0.0) {
      /* 2 T_n(x)/n */
      *result = 2.0 * cos(n * acos(x)) / n;
      return GSL_SUCCESS;
    }
    else {
      int k;
      double gkm2 = gsl_sf_gegenpoly_2(lambda, x);
      double gkm1 = gsl_sf_gegenpoly_3(lambda, x);
      double gk;
      for(k=4; k<=n; k++) {
        gk = (2.0*(k+lambda-1.0)*x*gkm1 - (k+2.0*lambda-2.0)*gkm2) / k;
	gkm2 = gkm1;
	gkm1 = gk;
      }
      *result = gk;
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_gegenpoly_array_impl(int nmax, double lambda, double x, double * result_array)
{
  int k;

  if(lambda <= -0.5 || nmax < 0) return GSL_EDOM;

  /* n == 0 */
  result_array[0] = 1.0;
  if(nmax == 0) return GSL_SUCCESS;

  /* n == 1 */
  if(lambda == 0.0)
    result_array[1] = 2.0*x;
  else
    result_array[1] = 2.0*lambda*x;

  /* n <= nmax */
  for(k=2; k<=nmax; k++) {
    double term1 = 2.0*(k+lambda-1.0) * x * result_array[k-1];
    double term2 = (k+2.0*lambda-2.0)	  * result_array[k-2];
    result_array[k] = (term1 - term2) / k;
  }
  
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_gegenpoly_n_e(int n, double lambda, double x, double * result)
{
  int status = gsl_sf_gegenpoly_n_impl(n, lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gegenpoly_n_e", status);
  }
  return status;
}

int
gsl_sf_gegenpoly_array_e(int nmax, double lambda, double x, double * result_array)
{
  int status = gsl_sf_gegenpoly_n_impl(nmax, lambda, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gegenpoly_array_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_gegenpoly_n(int n, double lambda, double x)
{
  double y;
  int status = gsl_sf_gegenpoly_n_impl(n, lambda, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_gegenpoly_n", status);
  }
  return y;
}
