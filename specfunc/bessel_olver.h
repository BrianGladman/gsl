/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef BESSEL_OLVER_H_
#define BESSEL_OLVER_H_


int gsl_sf_bessel_Jnu_asymp_Olver_impl(double nu, double x, double * result);
int gsl_sf_bessel_Ynu_asymp_Olver_impl(double nu, double x, double * result);


#endif  /* !BESSEL_OLVER_H_ */
