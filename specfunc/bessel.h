/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef _BESSEL_H_
#define _BESSEL_H_

int gsl_sf_bessel_Inu_Jnu_taylor_impl(double nu, double x,
                                      int sign,
                                      int kmax,
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

int gsl_sf_bessel_Jnu_asymp_trans_impl(double nu, double x, double * result);
int gsl_sf_bessel_Ynu_asymp_trans_impl(double nu, double x, double * result);

int gsl_sf_bessel_Jnu_asymp_Olver_impl(double nu, double x, double * result);
int gsl_sf_bessel_Ynu_asymp_Olver_impl(double nu, double x, double * result);


#endif /* !_BESSEL_H_ */
