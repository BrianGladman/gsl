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


int
gsl_sf_bessel_JnuYnu_zero(double nu, double x,
                          double * Jnu,  double * Ynu,
                          double * Jpnu, double * Ypnu
		          );


int
gsl_sf_bessel_J_recur(double nu_min, double x, const int kmax,
                      double J_start, double Jp_start,
	              double * J_end, double * Jp_end,
	              double * J_array, double * Jp_array);

int
gsl_sf_bessel_Y_recur(double nu_min, double x, const int kmax,
                      double Y_start, double Yp_start,
		      double * Y_end, double * Yp_end,
                      double * Y_array, double * Yp_array);

int
gsl_sf_bessel_I_recur(double nu_min, double x, const int kmax,
                      double I_start, double Ip_start,
	              double * I_end, double * Ip_end,
	              double * I_array, double * Ip_array);

int
gsl_sf_bessel_K_recur(double nu_min, double x, const int kmax,
                      double K_start, double Kp_start,
		      double * K_end, double * Kp_end,
                      double * K_array, double * Kp_array);

int
gsl_sf_bessel_JY_steed_CF2(double nu, double x,
                           double * P, double * Q);


#endif /* !_BESSEL_H_ */
