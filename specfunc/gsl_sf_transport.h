/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_TRANSPORT_H__
#define __GSL_SF_TRANSPORT_H__

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


/* Transport function:
 *   J(n,x) := Integral[ t^n e^t /(e^t - 1)^2, {t,0,x}]
 */

/* J(2,x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_transport_2_impl(const double x, gsl_sf_result * result);
int     gsl_sf_transport_2_e(const double x, gsl_sf_result * result);


/* J(3,x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_transport_3_impl(const double x, gsl_sf_result * result);
int     gsl_sf_transport_3_e(const double x, gsl_sf_result * result);


/* J(4,x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_transport_4_impl(const double x, gsl_sf_result * result);
int     gsl_sf_transport_4_e(const double x, gsl_sf_result * result);


/* J(5,x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_transport_5_impl(const double x, gsl_sf_result * result);
int     gsl_sf_transport_5_e(const double x, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_TRANSPORT_H__ */
