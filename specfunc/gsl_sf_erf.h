/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_ERF_H_
#define GSL_SF_ERF_H_

/* erfc(x); safe for all x */
double gsl_sf_erfc(double x);
double gsl_sf_log_erfc(double x);

/* erf(x); safe for all x */
double gsl_sf_erf(double x);

/* Probability functions:
 * Z(x) :  Abramowitz+Stegun 26.2.1
 * Q(x) :  Abramowitz+Stegun 26.2.3
 * safe for all x
 */
double gsl_sf_erf_Z(double x);
double gsl_sf_erf_Q(double x);

/* Asymptotic expansion of erfc(x) for x -> +infinity.
 * These are not safe for all x.
 */
double gsl_sf_erfc_asymptotic(double x);
double gsl_sf_log_erfc_asymptotic(double x);


#endif /* !GSL_SF_ERF_H_ */
