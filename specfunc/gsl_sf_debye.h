/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_DEBYE_H__
#define __GSL_SF_DEBYE_H__

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


/* D_n(x) := n/x^n Integrate[t^n/(e^t - 1), {t,0,x}] */

/* D_1(x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_debye_1_impl(const double x, gsl_sf_result * result);
int     gsl_sf_debye_1_e(const double x, gsl_sf_result * result);


/* D_2(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_2_impl(const double x, gsl_sf_result * result);
int     gsl_sf_debye_2_e(const double x, gsl_sf_result * result);


/* D_3(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_3_impl(const double x, gsl_sf_result * result);
int     gsl_sf_debye_3_e(const double x, gsl_sf_result * result);


/* D_4(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_4_impl(const double x, gsl_sf_result * result);
int     gsl_sf_debye_4_e(const double x, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_DEBYE_H__ */
