/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_COULOMB_H_
#define GSL_COULOMB_H_

#include <gsl_mode.h>
#include <gsl_sf_result.h>


/* Normalized hydrogenic bound states, radial dependence. */

/* R_1 := 2Z sqrt(Z) exp(-Z r)
 */
int gsl_sf_hydrogenicR_1_impl(double Z, double r, gsl_sf_result * result);
int gsl_sf_hydrogenicR_1_e(double Z, double r, gsl_sf_result * result);

/* R_n := norm exp(-Z r/n) (2Z/n)^l Laguerre[n-l-1, 2l+1, 2Z/n r]
 *
 * normalization such that psi(n,l,r) = R_n Y_{lm}
 */
int gsl_sf_hydrogenicR_impl(int n, int l, double Z, double r, gsl_sf_result * result);
int gsl_sf_hydrogenicR_e(int n, int l, double Z, double r, gsl_sf_result * result);


/* Coulomb wave functions F_{lam_F}(eta,x), G_{lam_G}(eta,x)
 * and their derivatives; lam_G := lam_F - k_lam_G
 *
 * lam_F, lam_G > -0.5
 * x > 0.0
 *
 * Conventions of Abramowitz+Stegun.
 *
 * Because their can be a large dynamic range of values,
 * overflows are handled gracefully. If an overflow occurs,
 * GSL_EOVRFLW is signalled and exponent(s) are returned
 * through exp_F, exp_G. These are such that
 *
 *   F_L(eta,x)  =  fc[k_L] * exp(exp_F)
 *   G_L(eta,x)  =  gc[k_L] * exp(exp_G)
 *   F_L'(eta,x) = fcp[k_L] * exp(exp_F)
 *   G_L'(eta,x) = gcp[k_L] * exp(exp_G)
 */
int
gsl_sf_coulomb_wave_FG_impl(const double eta, const double x,
                            const double lam_F,
			    const int  k_lam_G,
			    gsl_mode_t mode,
                            gsl_sf_result * F, gsl_sf_result * Fp,
			    gsl_sf_result * G, gsl_sf_result * Gp,
			    double * exp_F, double * exp_G);


/* F_L(eta,x)
 */
int gsl_sf_coulomb_wave_F_array_impl(double lam_min, int kmax,
                               double eta, double x,
			       gsl_mode_t mode,
                               gsl_sf_result * fc,
			       double * F_exponent
                               );
int gsl_sf_coulomb_wave_F_e(double lam_min, int kmax,
                            double eta, double x,
			    gsl_mode_t mode,
                            gsl_sf_result * fc,
			    double * F_exponent
                            );

/* F_L(eta,x), G_L(eta,x)
 */
int gsl_sf_coulomb_wave_FG_array_impl(double lam_min, int kmax,
                                double eta, double x,
			        gsl_mode_t mode,
                                double * fc_array, double * gc_array,
				double * F_exponent,
				double * G_exponent
                                );
int gsl_sf_coulomb_wave_FG_array_e(double lam_min, int kmax,
                             double eta, double x,
                             double * fc_array, double * gc_array,
			     double * F_exponent,
			     double * G_exponent,
			     gsl_mode_t mode, unsigned int err_bits
                             );

/* F_L(eta,x), G_L(eta,x), F'_L(eta,x), G'_L(eta,x)
 */
int gsl_sf_coulomb_wave_FGp_impl(double lam_min, int kmax,
                                 double eta, double x,
			         gsl_mode_t mode,
                                 gsl_sf_result * fc, gsl_sf_result * fcp,
                                 gsl_sf_result * gc, gsl_sf_result * gcp,
				 double * F_exponent,
				 double * G_exponent
	      	      	         );
int gsl_sf_coulomb_wave_FGp_e(double lam_min, int kmax,
                              double eta, double x,
			      gsl_mode_t mode,
                              gsl_sf_result * fc, gsl_sf_result * fcp,
                              gsl_sf_result * gc, gsl_sf_result * gcp,
			      double * F_exponent,
			      double * G_exponent
	      	      	      );

/* Coulomb wave function divided by the argument,
 * F(xi, eta)/xi. This is the function which reduces to
 * spherical Bessel functions in the limit eta->0.
 */
int gsl_sf_coulomb_wave_sphF_array_impl(double lam_min, int kmax,
                                        double eta, double x,
					gsl_mode_t mode,
	              	                double * fc_array,
			                double * F_exponent
	      	      	                );


/* Coulomb wave function normalization constant.
 * [Abramowitz+Stegun 14.1.8, 14.1.9]
 */
int gsl_sf_coulomb_CL_impl(double L, double eta, double * result);
int gsl_sf_coulomb_CL_list(double Lmin, int kmax, double eta, double *cl);


#endif  /* !GSL_COULOMB_H_ */
