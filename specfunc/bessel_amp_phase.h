/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef _BESSEL_AMP_PHASE_H_
#define _BESSEL_AMP_PHASE_H_

#include "gsl_sf_chebyshev.h"

extern const gsl_sf_cheb_series _bessel_amp_phase_bm0_cs;
extern const gsl_sf_cheb_series _bessel_amp_phase_bth0_cs;

extern const gsl_sf_cheb_series _bessel_amp_phase_bm1_cs;
extern const gsl_sf_cheb_series _bessel_amp_phase_bth1_cs;


/* large argument expansions [Abramowitz+Stegun, 9.2.28-29]; x > 0 */
double gsl_sf_bessel_asymp_Mnu(double nu, double x);
double gsl_sf_bessel_asymp_thetanu(double nu, double x);


#endif /* !_BESSEL_AMP_PHASE_H_ */
