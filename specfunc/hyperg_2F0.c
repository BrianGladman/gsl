/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "hyperg.h"
#include "gsl_sf_hyperg.h"


int
gsl_sf_hyperg_2F0_impl(const double a, const double b, const double x, gsl_sf_result * result)
{
  if(x < 0.0) {
    /* Use "definition" 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x).
     */
    gsl_sf_result U;
    double pre = pow(-1.0/x, a);
    int stat_U = gsl_sf_hyperg_U_impl(a, 1.0+a-b, -1.0/x, &U);
    result->val = pre * U.val;
    result->err = GSL_DBL_EPSILON * fabs(result->val) + pre * U.err;
    return stat_U;
  }
  else if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* Use asymptotic series. ??
     */
    /* return hyperg_2F0_series(a, b, x, -1, result, &prec); */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
}


int
gsl_sf_hyperg_2F0_e(const double a, const double b, const double x, gsl_sf_result * result)
{
  int status = gsl_sf_hyperg_2F0_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_2F0_e", status);
  }
  return status;
}
