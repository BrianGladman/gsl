/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_LAGUERRE_H_
#define GSL_SF_LAGUERRE_H_

#include <gsl_sf_result.h>


/* L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x) */


/* Evaluate generalized Laguerre polynomials
 * using explicit representations.
 *
 * exceptions: none
 */
int gsl_sf_laguerre_1_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_laguerre_2_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_laguerre_3_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_laguerre_1_e(double a, double x, gsl_sf_result * result);
int gsl_sf_laguerre_2_e(double a, double x, gsl_sf_result * result);
int gsl_sf_laguerre_3_e(double a, double x, gsl_sf_result * result);


/* Evaluate generalized Laguerre polynomials.
 *
 * a > -1.0
 * n >= 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_laguerre_n_impl(const int n, const double a, const double x, gsl_sf_result * result);
int     gsl_sf_laguerre_n_e(int n, double a, double x, gsl_sf_result * result);


#endif /* GSL_SF_LAGUERRE_H_ */
