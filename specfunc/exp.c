/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"


int gsl_sf_exp_impl(double x, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = exp(x);
    return GSL_SUCCESS;
  }
}
