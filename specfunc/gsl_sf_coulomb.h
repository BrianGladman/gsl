/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_COULOMB_H_
#define GSL_COULOMB_H_


/* Normalized hydrogenic bound states, radial dependence.
 */
int    gsl_sf_hydrogenicR_1_impl(double Z, double r, double * result);
int    gsl_sf_hydrogenicR_2_impl(int l, double Z, double r, double * result);
int    gsl_sf_hydrogenicR_impl(int n, int l, double Z, double r, double * result);
int    gsl_sf_hydrogenicR_1_e(double Z, double r, double * result);
int    gsl_sf_hydrogenicR_2_e(int l, double Z, double r, double * result);
int    gsl_sf_hydrogenicR_e(int n, int l, double Z, double r, double * result);
double gsl_sf_hydrogenicR_1(double Z, double r);
double gsl_sf_hydrogenicR_2(int l, double Z, double r);
double gsl_sf_hydrogenicR(int n, int l, double Z, double r);



/* Coulomb wave functions F_L(eta,x), G_L(eta,x) and their
 * derivatives. Evaluates for a squence of (real-valued) L:
 *
 *    { L } = { lam_min, lam_min+1, ..., lam_min+kmax }
 *
 *    { F } = { fc[0], fc[1], ..., fc[kmax] }
 *    { G } = { gc[0], gc[1], ..., gc[kmax] }
 *     etc.
 *
 * lam_min > -0.5
 * kmax >= 0
 * x > 0.0
 *
 * Conventions of Abramowitz+Stegun.
 *
 * Because their can be a large dynamic range of values,
 * overflows are handled gracefully. If an overflow occurs,
 * GSL_EOVRFLW is signalled and exponent(s) are returned
 * through F_exponent, G_exponent. These are such that
 *
 *   F_L(eta,x)  =  fc[k_L] * exp(F_exponent)
 *   G_L(eta,x)  =  gc[k_L] * exp(G_exponent)
 *   F_L'(eta,x) = fcp[k_L] * exp(F_exponent)
 *   G_L'(eta,x) = gcp[k_L] * exp(G_exponent)
 */

int gsl_sf_coulomb_wave_F_impl(double lam_min, int kmax,
                               double eta, double x,
                               double * fc,
			       double * F_exponent
                               );
int gsl_sf_coulomb_wave_F_e(double lam_min, int kmax,
                            double eta, double x,
                            double * fc,
			    double * F_exponent
                            );

int gsl_sf_coulomb_wave_FG_impl(double lam_min, int kmax,
                                double eta, double x,
                                double * fc, double * gc,
				double * F_exponent,
				double * G_exponent
                                );
int gsl_sf_coulomb_wave_FG_e(double lam_min, int kmax,
                             double eta, double x,
                             double * fc, double * gc,
			     double * F_exponent,
			     double * G_exponent
                             );

int gsl_sf_coulomb_wave_FGp_impl(double lam_min, int kmax,
                                 double eta, double x,
                                 double * fc, double * fcp,
                                 double * gc, double * gcp,
				 double * F_exponent,
				 double * G_exponent
	      	      	         );
int gsl_sf_coulomb_wave_FGp_e(double lam_min, int kmax,
                              double eta, double x,
                              double * fc, double * fcp,
                              double * gc, double * gcp,
			      double * F_exponent,
			      double * G_exponent
	      	      	      );

/* Coulomb wave function divided by the argument,
 * F(xi, eta)/xi. This is the function which reduces to
 * spherical Bessel functions in the limit eta->0.
 */
int gsl_sf_coulomb_wave_sphF_impl(double lam_min, int kmax,
                                  double eta, double x,
	              	          double * fc,
			          double * F_exponent
	      	      	          );


/* Coulomb wave function normalization constant.
 * [Abramowitz+Stegun 14.1.8, 14.1.9]
 */
int    gsl_sf_coulomb_CL_impl(double L, double eta, double * result);
int    gsl_sf_coulomb_CL_list(double Lmin, int kmax, double eta, double *cl);


#endif  /* !GSL_COULOMB_H_ */
