/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_LOG_H__
#define __GSL_SF_LOG_H__

#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS


/* Provide a logarithm function with GSL semantics.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_impl(const double x, gsl_sf_result * result);
int gsl_sf_log_e(const double x, gsl_sf_result * result);


/* Log(|x|)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_abs_impl(const double x, gsl_sf_result * result);
int gsl_sf_log_abs_e(const double x, gsl_sf_result * result);


/* Complex Logarithm
 *   exp(lnr + I theta) = zr + I zi
 * Returns argument in [-pi,pi].
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_impl(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta);
int gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta);


/* Log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_impl(const double x, gsl_sf_result * result);
int gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result);


/* Log(1 + x) - x
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_mx_impl(const double x, gsl_sf_result * result);
int gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result);


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
#endif /* HAVE_INLINE */


__END_DECLS

#endif /* __GSL_SF_LOG_H__ */
