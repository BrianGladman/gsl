/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_EXP_H_
#define GSL_SF_EXP_H_


/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 */
int gsl_sf_exp_impl(double x, double * result);


#ifdef HAVE_INLINE
#include <gsl_errno.h>
extern inline
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
  else
    *result = exp(x);
    return GSL_SUCCESS;
  }  
}
#endif  /* HAVE_INLINE */


#endif  /* !GSL_SF_EXP_H_ */
