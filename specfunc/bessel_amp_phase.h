/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef _BESSEL_AMP_PHASE_H_
#define _BESSEL_AMP_PHASE_H_

#include "gsl_sf_chebyshev.h"

const struct gsl_sf_ChebSeries 	_bessel_amp_phase_bm0_cs;
const struct gsl_sf_ChebSeries   _bessel_amp_phase_bth0_cs;

const struct gsl_sf_ChebSeries   _bessel_amp_phase_bm1_cs;
const struct gsl_sf_ChebSeries   _bessel_amp_phase_bth1_cs;


/* large argument expansions [Abramowitz+Stegun, 9.2.28-29]; x > 0 */

inline double gsl_sf_bessel_asymp_Mnu(double nu, double x)
{
  double x_inv  = 1./x;
  double x_inv2 = x_inv*x_inv;
  double mu     = 4.*nu*nu;
  double Mnu2 = 2./(M_PI) * x_inv * (1. + (mu-1.)/8.*x_inv2
                                        + (mu-1.)*(mu-9.)*3./128.*x_inv2*x_inv2
                                    );
  return sqrt(Mnu2);
}

inline double gsl_sf_bessel_asymp_thetanu(double nu, double x)
{
  double x_inv  = 1./x;
  double x_inv2 = x_inv*x_inv;
  double mu     = 4.*nu*nu;
  return x * (1. - (0.5*nu + 0.25)*M_PI*x_inv 
	         + (mu-1.)/8. * x_inv2
                 + (mu-1.)*(mu-25.)/384. *x_inv2*x_inv2
	     );
}


#endif /* !_BESSEL_AMP_PHASE_H_ */
