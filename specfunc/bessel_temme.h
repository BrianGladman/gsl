/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef BESSEL_TEMME_H_
#define BESSEL_TEMME_H_

#include "gsl_sf_result.h"


int
gsl_sf_bessel_Y_temme(double nu, double x,
                      gsl_sf_result * Y_nu,
                      gsl_sf_result * Y_nup1);

int
gsl_sf_bessel_K_scaled_temme(double nu, double x,
                             double * K_nu, double * K_nup1, double * Kp_nu);


#endif  /* !BESSEL_TEMME_H_ */
