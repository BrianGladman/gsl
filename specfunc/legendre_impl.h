/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef LEGENDRE_IMPL_H_
#define LEGENDRE_IMPL_H_


int gsl_sf_legendre_Pl_impl(int l, double x, double * result, double * harvest);

int gsl_sf_legendre_Plm_impl(int l, int m, double one_m_x, double one_p_x, double * result, double * harvest);

int gsl_sf_legendre_Plm_impl(int l, int m, double one_m_x, double one_p_x, double * result, double * harvest);

int gsl_sf_conical_xlt1_large_mu_impl(double mu, double tau, double x, double * result);


#endif /* !LEGENDRE_IMPL_H_ */
