/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_ERF_H_
#define GSL_SF_ERF_H_


/* Complementary Error Function
 * erfc(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,x,Infinity}]
 *
 * exceptions: none
 */
double gsl_sf_erfc(double x);
double gsl_sf_log_erfc(double x);


/* Error Function
 * erf(x) := 2/Sqrt[Pi] Integrate[Exp[-t^2], {t,0,x}]
 *
 * exceptions: none
 */
double gsl_sf_erf(double x);


/* Probability functions:
 * Z(x) :  Abramowitz+Stegun 26.2.1
 * Q(x) :  Abramowitz+Stegun 26.2.3
 *
 * exceptions: none
 */
double gsl_sf_erf_Z(double x);
double gsl_sf_erf_Q(double x);


/* Asymptotic expansion of erfc(x) for x -> +infinity.
 *
 */
double gsl_sf_erfc_asymptotic(double x);
double gsl_sf_log_erfc_asymptotic(double x);


#endif /* !GSL_SF_ERF_H_ */
