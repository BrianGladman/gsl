/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_DAWSON_H__
#define __GSL_SF_DAWSON_H__

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


/* Dawson's integral:
 *
 *   Exp[-x^2] Integral[ Exp[t^2], {t,0,x}]
 *
 * exceptions: GSL_EUNDRFLW;
 */
int     gsl_sf_dawson_impl(double x, gsl_sf_result * result);
int     gsl_sf_dawson_e(double x, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_DAWSON_H__ */
