/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Miscellaneous implementations of use
 * for evaluation of hypergeometric functions.
 */
#ifndef _HYPERG_H_
#define _HYPERG_H_


/* Direct implementation of 1F1 series.
 */
int
gsl_sf_hyperg_1F1_series_impl(double a, double b, double x, double * result, double * prec);


/* Implementation of the 1F1 related to the
 * incomplete gamma function: 1F1(1,b,x), b >= 1.
 */
int
gsl_sf_hyperg_1F1_1_impl(double b, double x, double * result);


/* 1F1(1,b,x) for integer b >= 1
 */
int
gsl_sf_hyperg_1F1_1_int_impl(int b, double x, double * result);


/* Implementation of large b asymptotic.
 * [Bateman v. I, 6.13.3 (18)]
 * [Luke, The Special Functions and Their Approximations v. I, p. 129, 4.8 (4)]
 *
 * a^2 << b, |x|/|b| < 1 - delta
 */
int
gsl_sf_hyperg_1F1_large_b_impl(double a, double b, double x, double * result);


/* Implementation of large b asymptotic.
 *
 * Assumes a > 0 is small, x > 0, and |x|<|b|.
 */
int
gsl_sf_hyperg_U_large_b_impl(double a, double b, double x,
                             double * result,
                             double * ln_multiplier
                             );


/* Implementation of 2F0 asymptotic series.
 */
int
gsl_sf_hyperg_2F0_series_impl(double a, double b, double x, int n_trunc,
                              double * result, double * prec);


#endif  /* !_HYPERG_H_ */
