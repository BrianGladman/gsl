/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_GEGENBAUER_H_
#define GSL_SF_GEGENBAUER_H_


/* Evaluate Gegenbauer polynomials
 * using explicit representations.
 */
double gsl_sf_gegenpoly_1(double lambda, double x);
double gsl_sf_gegenpoly_2(double lambda, double x);
double gsl_sf_gegenpoly_3(double lambda, double x);


/* Evaluate Gegenbauer polynomials
 * by recurrence on n.
 */

int gsl_sf_gegenpoly_n_e(int n, double lambda, double x, double * result);

double gsl_sf_gegenpoly_n(int n, double lambda, double x);

int gsl_sf_gegenpoly_n_impl(int n, double lambda, double x, double * result);


/* Calculate array of Gegenbauer polynomials
 * for a n = (0, 1, 2, ... nmax)
 */



int gsl_sf_gegenpoly_array_impl(int nmax, double lambda, double x, double * result_array);


#endif  /* !GSL_SF_GEGENBAUER_H_ */
