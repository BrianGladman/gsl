/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_SYNCHROTRON_H__
#define __GSL_SF_SYNCHROTRON_H__

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


/* First synchrotron function:
 *   synchrotron_1(x) = x Integral[ K_{5/3}(t), {t, x, Infinity}]
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_synchrotron_1_impl(double x, gsl_sf_result * result);
int     gsl_sf_synchrotron_1_e(double x, gsl_sf_result * result);


/* Second synchroton function:
 *   synchrotron_2(x) = x * K_{2/3}(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_synchrotron_2_impl(double x, gsl_sf_result * result);
int     gsl_sf_synchrotron_2_e(double x, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_SYNCHROTRON_H__ */
