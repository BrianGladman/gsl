/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_ERF_H__
#define __GSL_SF_ERF_H__

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


/* Complementary Error Function
 * erfc(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,x,Infinity}]
 *
 * exceptions: none
 */
int gsl_sf_erfc_impl(double x, gsl_sf_result * result);
int gsl_sf_erfc_e(double x, gsl_sf_result * result);


/* Log Complementary Error Function
 *
 * exceptions: none
 */
int gsl_sf_log_erfc_impl(double x, gsl_sf_result * result);
int gsl_sf_log_erfc_e(double x, gsl_sf_result * result);


/* Error Function
 * erf(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,0,x}]
 *
 * exceptions: none
 */
int gsl_sf_erf_impl(double x, gsl_sf_result * result);
int gsl_sf_erf_e(double x, gsl_sf_result * result);


/* Probability functions:
 * Z(x) :  Abramowitz+Stegun 26.2.1
 * Q(x) :  Abramowitz+Stegun 26.2.3
 *
 * exceptions: none
 */
int gsl_sf_erf_Z_impl(double x, gsl_sf_result * result);
int gsl_sf_erf_Q_impl(double x, gsl_sf_result * result);
int gsl_sf_erf_Z_e(double x, gsl_sf_result * result);
int gsl_sf_erf_Q_e(double x, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_ERF_H__ */
