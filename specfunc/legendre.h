/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Declare private but non-local support functions
 * used in various Legendre function evaluations.
 */

#include <gsl/gsl_sf_result.h>


/* Large negative mu asymptotic
 * P^{-mu}_{-1/2 + I tau}, mu -> Inf
 * |x| < 1
 */
int
gsl_sf_conicalP_xlt1_large_neg_mu_impl(double mu, double tau, double x,
                                       gsl_sf_result * result, double * ln_multiplier);


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}, tau -> Inf
 * 1 < x
 */
int
gsl_sf_conicalP_xgt1_neg_mu_largetau_impl(const double mu, const double tau,
                                          const double x, double acosh_x,
                                          gsl_sf_result * result, double * ln_multiplier);


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}, tau -> Inf 
 * -1 < x < 1
 */
int
gsl_sf_conicalP_xlt1_neg_mu_largetau_impl(const double mu, const double tau,
                                          const double x, const double acos_x,
                                          gsl_sf_result * result, double * ln_multiplier);


/* P^{mu}_{-1/2 + I tau}
 * x->Inf
 *
 *  * This is effective to precision EPS for
 *
 *    (mu^2 + tau^2)/((1 + tau^2)^(1/2) x^2) < EPS^{1/3}
 *
 * since it goes only to a fixed order, based on the
 * representation in terms of hypegeometric functions
 * of argument 1/x^2.
 * [Zhurina+Karmazina, (3.8)]
 */
int
gsl_sf_conicalP_large_x_impl(const double mu, const double tau, const double x,
                             gsl_sf_result * result, double * ln_multiplier);
