/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gegenbauer.h"


/* See: [Thompson, Atlas for Computing Mathematical Functions] */


int
gsl_sf_gegenpoly_1_impl(double lambda, double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(lambda == 0.0) {
    result->val = 2.0*x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 2.0*lambda*x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_gegenpoly_2_impl(double lambda, double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(lambda == 0.0) {
    const double txx = 2.0*x*x;
    result->val = -1.0 + txx;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(txx);
    return GSL_SUCCESS;
  }
  else {
    result->val = lambda*(-1.0 + 2.0*(1.0+lambda)*x*x);
    result->err = GSL_DBL_EPSILON * (fabs(result->val) + fabs(lambda));
    return GSL_SUCCESS;
  }
}

int
gsl_sf_gegenpoly_3_impl(double lambda, double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(lambda == 0.0) {
    result->val = x*(-2.0 + 4.0/3.0*x*x);
    result->err = GSL_DBL_EPSILON * (fabs(result->val) + fabs(x));
    return GSL_SUCCESS;
  }
  else {
    double c = 4.0 + lambda*(6.0 + 2.0*lambda);
    result->val = 2.0*lambda * x * ( -1.0 - lambda + c*x*x/3.0 );
    result->err = GSL_DBL_EPSILON * (fabs(result->val) + fabs(lambda * x));
    return GSL_SUCCESS;
  }
}


int
gsl_sf_gegenpoly_n_impl(int n, double lambda, double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(lambda <= -0.5 || n < 0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(n == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    return gsl_sf_gegenpoly_1_impl(lambda, x, result);
  }
  else if(n == 2) {
    return gsl_sf_gegenpoly_2_impl(lambda, x, result);
  }
  else if(n == 3) {
    return gsl_sf_gegenpoly_3_impl(lambda, x, result);
  }
  else {
    if(lambda == 0.0 && (x >= -1.0 || x <= 1.0)) {
      /* 2 T_n(x)/n */
      const double z = n * acos(x);
      result->val = 2.0 * cos(z) / n;
      result->err = GSL_DBL_EPSILON * fabs(z * result->val);
      return GSL_SUCCESS;
    }
    else {
      int k;
      gsl_sf_result g2;
      gsl_sf_result g3;
      int stat_g2 = gsl_sf_gegenpoly_2_impl(lambda, x, &g2);
      int stat_g3 = gsl_sf_gegenpoly_3_impl(lambda, x, &g3);
      int stat_g  = GSL_ERROR_SELECT_2(stat_g2, stat_g3);
      double gkm2 = g2.val;
      double gkm1 = g3.val;
      double gk;
      for(k=4; k<=n; k++) {
        gk = (2.0*(k+lambda-1.0)*x*gkm1 - (k+2.0*lambda-2.0)*gkm2) / k;
	gkm2 = gkm1;
	gkm1 = gk;
      }
      result->val = gk;
      result->err = GSL_DBL_EPSILON * 0.5 * n * fabs(gk);
      return stat_g;
    }
  }
}


int
gsl_sf_gegenpoly_array_impl(int nmax, double lambda, double x, double * result_array)
{
  int k;

  if(result_array == 0) {
    return GSL_EFAULT;
  }

  if(lambda <= -0.5 || nmax < 0) {
    return GSL_EDOM;
  }

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
gsl_sf_gegenpoly_1_e(double lambda, double x, gsl_sf_result * result)
{
  int status = gsl_sf_gegenpoly_1_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gegenpoly_1_e", status);
  }
  return status;
}

int
gsl_sf_gegenpoly_2_e(double lambda, double x, gsl_sf_result * result)
{
  int status = gsl_sf_gegenpoly_2_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gegenpoly_2_e", status);
  }
  return status;
}

int
gsl_sf_gegenpoly_3_e(double lambda, double x, gsl_sf_result * result)
{
  int status = gsl_sf_gegenpoly_3_impl(lambda, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gegenpoly_3_e", status);
  }
  return status;
}

int
gsl_sf_gegenpoly_n_e(int n, double lambda, double x, gsl_sf_result * result)
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
  int status = gsl_sf_gegenpoly_array_impl(nmax, lambda, x, result_array);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_gegenpoly_array_e", status);
  }
  return status;
}
