/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_elementary.h"


int
gsl_sf_multiply_impl(const double x, const double y, gsl_sf_result * result)
{
  const double ax = fabs(x);
  const double ay = fabs(y);

  if(x == 0.0 || y == 0.0) {
    /* It is necessary to eliminate this immediately.
     */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if((ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0)) {
    /* Straddling 1.0 is always safe.
     */
    result->val = x*y;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double f = 1.0 - 2.0 * GSL_DBL_EPSILON;
    const double min = GSL_MIN_DBL(fabs(x), fabs(y));
    const double max = GSL_MAX_DBL(fabs(x), fabs(y));
    if(max < 0.9 * GSL_SQRT_DBL_MAX || min < (f * DBL_MAX)/max) {
      result->val = x*y;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return (fabs(result->val) == 0.0 ? GSL_EUNDRFLW : GSL_SUCCESS);
    }
    else {
      result->val = 0.0; /* FIXME: should be Inf */
      result->err = 0.0;
      return GSL_EOVRFLW;
    }
  }
}


int
gsl_sf_multiply_err_impl(const double x, const double dx,
                         const double y, const double dy,
                         gsl_sf_result * result)
{
  int status = gsl_sf_multiply_impl(x, y, result);
  result->err += fabs(dx*y) + fabs(dy*x);
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result)
{
  int status = gsl_sf_multiply_impl(x, y, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_multiply_e", status);
  }
  return status;
}


int
gsl_sf_multiply_err_e(const double x, const double dx,
                      const double y, const double dy,
                      gsl_sf_result * result)
{
  int status = gsl_sf_multiply_err_impl(x, dx, y, dy, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_multiply_err_e", status);
  }
  return status;
}
