/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef BESSEL_JN_IMPL_H_
#define BESSEL_JN_IMPL_H_

int gsl_sf_bessel_Jn_impl(int n, double x, double * result);
int gsl_sf_bessel_Jn_array_impl(int nmax, double x, double * result_array);

#endif /* ! BESSEL_JN_IMPL_H_ */
