/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_DILOG_H__
#define __GSL_SF_DILOG_H__

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


/* Real part of DiLogarithm(x), for real argument.
 * In Lewin's notation, this is Li_2(x).
 *
 *   Li_2(x) = - Re[ Integrate[ Log[1-s] / s, {s, 0, x}] ]
 *
 * Note that Im[Li_2(x)] = { 0 for x <= 1, -Pi*log(x) for x > 1 }
 */
int     gsl_sf_dilog_impl(const double x, gsl_sf_result * result);
int     gsl_sf_dilog_e(const double x, gsl_sf_result * result);


/* DiLogarithm(z), for complex argument z = r Exp[i theta].
 */
int gsl_sf_complex_dilog_impl(const double r, double theta, gsl_sf_result * result_re, gsl_sf_result * result_im);
int gsl_sf_complex_dilog_e(const double r, const double theta, gsl_sf_result * result_re, gsl_sf_result * result_im);


__END_DECLS

#endif /* __GSL_SF_DILOG_H__ */
