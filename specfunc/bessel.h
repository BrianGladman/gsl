/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef _BESSEL_H_
#define _BESSEL_H_

#include "gsl_sf_result.h"


int gsl_sf_bessel_Inu_Jnu_taylor_impl(double nu, double x,
                                      int sign,
                                      int kmax,
				      double threshold,
                                      double * result
                                      );

int gsl_sf_bessel_Jnu_asympx_impl(double nu, double x, gsl_sf_result * result);
int gsl_sf_bessel_Ynu_asympx_impl(double nu, double x, gsl_sf_result * result);

int gsl_sf_bessel_Inu_scaled_asympx_impl(double nu, double x, gsl_sf_result * result);
int gsl_sf_bessel_Knu_scaled_asympx_impl(double nu, double x, gsl_sf_result * result);

int gsl_sf_bessel_Inu_scaled_asymp_unif_impl(double nu, double x, gsl_sf_result * result);
int gsl_sf_bessel_Knu_scaled_asymp_unif_impl(double nu, double x, gsl_sf_result * result);


int
gsl_sf_bessel_JnuYnu_zero(double nu,
                          double * Jnu,  double * Ynu,
                          double * Jpnu, double * Ypnu
		          );

int
gsl_sf_bessel_J_CF1_ser(double nu, double x, double * ratio);

int
gsl_sf_bessel_I_CF1_ser(double nu, double x, double * ratio);


int
gsl_sf_bessel_JY_steed_CF2(double nu, double x,
                           double * P, double * Q);

int
gsl_sf_bessel_K_scaled_steed_temme_CF2(const double nu, const double x,
                                       double * K_nu, double * K_nup1,
				       double * Kp_nu);


/* These are of use in calculating the oscillating
 * Bessel functions.
 *   cos(y - pi/4 + eps)
 *   sin(y - pi/4 + eps)
 */
int gsl_sf_bessel_cos_pi4_impl(double y, double eps, gsl_sf_result * result);
int gsl_sf_bessel_sin_pi4_impl(double y, double eps, gsl_sf_result * result);


#endif /* !_BESSEL_H_ */
