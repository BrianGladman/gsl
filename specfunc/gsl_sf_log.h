/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_LOG_H_
#define GSL_SF_LOG_H_

#include <gsl/gsl_sf_result.h>


/* Provide a logarithm function with GSL semantics.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_impl(double x, gsl_sf_result * result);
int gsl_sf_log_e(double x, gsl_sf_result * result);


/* Log(|x|)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_abs_impl(double x, gsl_sf_result * result);
int gsl_sf_log_abs_e(double x, gsl_sf_result * result);


/* Complex Logarithm
 *   exp(lnr + I theta) = zr + I zi
 * Returns argument in [-pi,pi].
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_impl(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * theta);
int gsl_sf_complex_log_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * theta);


/* Log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_impl(double x, gsl_sf_result * result);
int gsl_sf_log_1plusx_e(double x, gsl_sf_result * result);


/* Log(1 + x) - x
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_mx_impl(double x, gsl_sf_result * result);
int gsl_sf_log_1plusx_mx_e(double x, gsl_sf_result * result);


#ifdef HAVE_INLINE
extern inline
int
gsl_sf_log_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x <= 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    result->val = log(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}
extern inline
int
gsl_sf_log_abs_impl(const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else {
    result->val = log(fabs(x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}
#endif  /* HAVE_INLINE */


#endif /* !GSL_SF_LOG_H_ */
