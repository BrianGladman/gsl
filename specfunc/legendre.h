/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Declare private but non-local support functions
 * used in various Legendre function evaluations.
 */


/* Large negative mu asymptotic
 * P^{-mu}_{-1/2 + I tau}, mu -> Inf
 * |x| < 1
 */
int
gsl_sf_conicalP_xlt1_large_neg_mu_impl(double mu, double tau, double x, double * result);


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}, tau -> Inf
 * 1 < x
 */
int
gsl_sf_conicalP_xgt1_neg_mu_largetau_impl(const double mu, const double tau,
                                          const double x, double * result);


/* Large tau uniform asymptotics
 * P^{-mu}_{-1/2 + I tau}, tau -> Inf 
 * -1 < x < 1
 */
int
gsl_sf_conicalP_xlt1_neg_mu_largetau_impl(const double mu, const double tau,
                                          const double x, double * result);

