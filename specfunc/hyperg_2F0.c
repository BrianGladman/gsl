/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <gsl_math.h>
#include <gsl_errno.h>
#include "hyperg.h"
#include "gsl_sf_hyperg.h"


int
gsl_sf_hyperg_2F0_impl(const double a, const double b, const double x, double * result)
{
  if(x < 0.0) {
    /* Use "definition" 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x).
     */
    double U;
    double pre = pow(-1.0/x, a);
    int stat_U = gsl_sf_hyperg_U_impl(a, 1.0+a-b, -1.0/x, &U);
    *result = pre * U;
    return stat_U;
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
    /* Use asymptotic series. ??
     */
    /* return hyperg_2F0_series(a, b, x, -1, result, &prec); */
    return GSL_EDOM;
  }
}


int
gsl_sf_hyperg_2F0_e(const double a, const double b, const double x, double * result)
{
  int status = gsl_sf_hyperg_2F0_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_2F0_e", status);
  }
  return status;
}

double
gsl_sf_hyperg_2F0(const double a, const double b, const double x)
{
  double y;
  int status = gsl_sf_hyperg_2F0_impl(a, b, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_2F0", status);
  }
  return y;
}
