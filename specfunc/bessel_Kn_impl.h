/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef BESSEL_KN_IMPL_H_
#define BESSEL_KN_IMPL_H_

int gsl_sf_bessel_Kn_scaled_impl(int n, double x, double * result);
int gsl_sf_bessel_Kn_impl(int n, double x, double * result);

int gsl_sf_bessel_Kn_scaled_array_impl(int nmax, double x, double * result_array);
int gsl_sf_bessel_Kn_array_impl(int nmax, double x, double * result_array);

#endif /* ! BESSEL_KN_IMPL_H_ */
