/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_LAGUERRE_H_
#define GSL_SF_LAGUERRE_H_


/* Evaluate generalized Laguerre polynomials
 * using explicit representations.
 *
 * exceptions: none
 */
double gsl_sf_laguerre_1(double a, double x);
double gsl_sf_laguerre_2(double a, double x);
double gsl_sf_laguerre_3(double a, double x);
double gsl_sf_laguerre_4(double a, double x);
double gsl_sf_laguerre_5(double a, double x);


/* Evaluate generalized Laguerre polynomials
 * using continued product representation.
 * For n<=5 it is better to use the explicit
 * functions, which are faster by approximately
 * a factor of 2--3.
 *
 * a > -1.0
 * n >= 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_laguerre_n_impl(const int n, const double a, const double x, double * result);
int     gsl_sf_laguerre_n_e(int n, double a, double x, double * result);
double  gsl_sf_laguerre_n(int n, double a, double x);


#endif /* GSL_SF_LAGUERRE_H_ */
