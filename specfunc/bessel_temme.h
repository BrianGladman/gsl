/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef BESSEL_TEMME_H_
#define BESSEL_TEMME_H_


int
gsl_sf_bessel_Y_temme(double nu, double x,
                      double * Y_nu, double * Y_nup1, double * Yp_nu);

int
gsl_sf_bessel_K_scaled_temme(double nu, double x,
                             double * K_nu, double * K_nup1, double * Kp_nu);


#endif  /* !BESSEL_TEMME_H_ */
