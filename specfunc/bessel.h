/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef _BESSEL_H_
#define _BESSEL_H_

int gsl_sf_bessel_Inu_Jnu_taylor_impl(double nu, double x,
                                      int sign,
                                      int kmax,
				      double threshold,
                                      double * result
                                      );

int gsl_sf_bessel_Jnu_asympx_impl(double nu, double x, double * result);
int gsl_sf_bessel_Ynu_asympx_impl(double nu, double x, double * result);

int gsl_sf_bessel_Inu_scaled_asympx_impl(double nu, double x, double * result);
int gsl_sf_bessel_Knu_scaled_asympx_impl(double nu, double x, double * result);

int gsl_sf_bessel_Inu_scaled_asymp_unif_impl(double nu, double x, double * result);
int gsl_sf_bessel_Knu_scaled_asymp_unif_impl(double nu, double x, double * result);

int gsl_sf_bessel_Jnu_asymp_Debye_impl(double nu, double x, double * result);
int gsl_sf_bessel_Ynu_asymp_Debye_impl(double nu, double x, double * result);

int gsl_sf_bessel_Jnu_asymp_Debye_osc_impl(double nu, double x, double * result);
int gsl_sf_bessel_Ynu_asymp_Debye_osc_impl(double nu, double x, double * result);


int
gsl_sf_bessel_JnuYnu_zero(double nu,
                          double * Jnu,  double * Ynu,
                          double * Jpnu, double * Ypnu
		          );

int
gsl_sf_bessel_JY_steed_CF2(double nu, double x,
                           double * P, double * Q);

int
gsl_sf_bessel_K_scaled_steed_temme_CF2(const double nu, const double x,
                                       double * K_nu, double * K_nup1,
				       double * Kp_nu);

#endif /* !_BESSEL_H_ */
