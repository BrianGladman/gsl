/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef BESSEL_TEMME_H_
#define BESSEL_TEMME_H_


int
gsl_sf_temme_gamma(double nu,
                   double * g_1pnu, double * g_1mnu,
                   double * g1, double * g2);


#endif  /* !BESSEL_TEMME_H_ */
