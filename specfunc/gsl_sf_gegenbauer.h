/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_GEGENBAUER_H
#define GSL_SF_GEGENBAUER_H

#include <gsl_sf_result.h>


/* Evaluate Gegenbauer polynomials
 * using explicit representations.
 *
 * exceptions: none
 */
int gsl_sf_gegenpoly_1_impl(double lambda, double x, gsl_sf_result * result);
int gsl_sf_gegenpoly_2_impl(double lambda, double x, gsl_sf_result * result);
int gsl_sf_gegenpoly_3_impl(double lambda, double x, gsl_sf_result * result);
int gsl_sf_gegenpoly_1_e(double lambda, double x, gsl_sf_result * result);
int gsl_sf_gegenpoly_2_e(double lambda, double x, gsl_sf_result * result);
int gsl_sf_gegenpoly_3_e(double lambda, double x, gsl_sf_result * result);


/* Evaluate Gegenbauer polynomials.
 *
 * lambda > -1/2, n >= 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_gegenpoly_n_impl(int n, double lambda, double x, gsl_sf_result * result);
int gsl_sf_gegenpoly_n_e(int n, double lambda, double x, gsl_sf_result * result);


/* Calculate array of Gegenbauer polynomials
 * for n = (0, 1, 2, ... nmax)
 *
 * lambda > -1/2, nmax >= 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_gegenpoly_array_impl(int nmax, double lambda, double x, double * result_array);
int gsl_sf_gegenpoly_array_e(int nmax, double lambda, double x, double * result_array);


#endif  /* !GSL_SF_GEGENBAUER_H */
