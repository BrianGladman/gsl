/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GAMMA_IMPL_H_
#define GAMMA_IMPL_H_

int gsl_sf_lngamma_impl(double x, double * result);
int gsl_sf_lngamma_complex_impl(double zr, double zi, double * lnr, double * arg);

int gsl_sf_fact_impl(int n, double * result);
int gsl_sf_lnfact_impl(int n, double * result);

int gsl_sf_choose_impl(unsigned int n, unsigned int m, double * result);
int gsl_sf_lnchoose_impl(unsigned int n, unsigned int m, double * result);


#endif /* !GAMMA_IMPL_H_ */
