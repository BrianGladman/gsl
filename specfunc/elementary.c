/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_elementary.h"


int
gsl_sf_multiply_impl(const double x, const double y, double * result)
{
  const double ax = fabs(x);
  const double ay = fabs(y);

  if(x == 0.0 || y == 0.0) {
    /* It is necessary to eliminate this immediately.
     */
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if((ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0)) {
    /* Straddling 1.0 is always safe.
     */
    *result = x*y;
    return GSL_SUCCESS;
  }
  else if(   (ax < 0.8*GSL_SQRT_DBL_MAX && ay < 0.8*GSL_SQRT_DBL_MAX)
          && (ax > 1.2*GSL_SQRT_DBL_MIN && ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    /* Not too big or too small. But just right...
     */
    *result = x*y;
    return GSL_SUCCESS;
  }
  else {
    const double lx = log(ax);
    const double ly = log(ay);
    if(lx+ly < GSL_LOG_DBL_MIN) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else if(lx+ly > GSL_LOG_DBL_MAX) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    else {
      /* Well... ok then.
       */
      *result = x*y;
      return GSL_SUCCESS;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_multiply_e(const double x, const double y, double * result)
{
  int status = gsl_sf_multiply_impl(x, y, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_multiply_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*/

double
gsl_sf_multiply(const double x, const double y)
{
  double r;
  int status = gsl_sf_multiply_impl(x, y, &r);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_multiply", status);
  }
  return r;
}
