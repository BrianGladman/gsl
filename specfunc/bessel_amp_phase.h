/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef _BESSEL_AMP_PHASE_H_
#define _BESSEL_AMP_PHASE_H_

#include "gsl_sf_chebyshev.h"

extern const gsl_sf_cheb_series _gsl_sf_bessel_amp_phase_bm0_cs;
extern const gsl_sf_cheb_series _gsl_sf_bessel_amp_phase_bth0_cs;

extern const gsl_sf_cheb_series _gsl_sf_bessel_amp_phase_bm1_cs;
extern const gsl_sf_cheb_series _gsl_sf_bessel_amp_phase_bth1_cs;


/* large argument expansions [Abramowitz+Stegun, 9.2.28-29]
 *
 * thetanu_corr = thetanu - x + 1/2 nu Pi
 *
 * assumes x > 0
 */
int gsl_sf_bessel_asymp_Mnu_impl(double nu, double x, double * result);
int gsl_sf_bessel_asymp_thetanu_corr_impl(double nu, double x, double * result); /* w/o x term */


#endif /* !_BESSEL_AMP_PHASE_H_ */
