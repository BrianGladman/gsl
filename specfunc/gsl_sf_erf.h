/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_ERF_H_
#define GSL_SF_ERF_H_

#include <gsl_sf_result.h>


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


#endif /* !GSL_SF_ERF_H_ */
