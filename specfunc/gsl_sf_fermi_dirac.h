/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_FERMI_DIRAC_H_
#define GSL_SF_FERMI_DIRAC_H_


/* Complete Fermi-Dirac Integrals:
 *
 *  F_j(x) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,0,Infinity}]
 *
 * Incomplete Fermi-Dirac Integrals:
 *
 *  F_j(x,b) := 1/Gamma[j+1] Integral[ t^j /(Exp[t-x] + 1), {t,b,Infinity}]
 */


/* Complete integral F_0(x) = ln(1 + e^x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_fermi_dirac_0_impl(double x, double * result);
int     gsl_sf_fermi_dirac_0_e(double x, double * result);
double  gsl_sf_fermi_dirac_0(double x);


/* Complete integral F_{-1}(x) = e^x / (1 + e^x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int     gsl_sf_fermi_dirac_m1_impl(double x, double * result);
int     gsl_sf_fermi_dirac_m1_e(double x, double * result);
double  gsl_sf_fermi_dirac_m1(double x);


/* Incomplete integral F_0(x,b) = ln(1 + e^(b-x)) - (b-x)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EDOM
 */
int     gsl_sf_fermi_dirac_inc_0_impl(double x, double b, double * result);
int     gsl_sf_fermi_dirac_inc_0_e(double x, double b, double * result);
double  gsl_sf_fermi_dirac_inc_0(double x, double b);



/* Fermi-Dirac function.
   zeta_inverse = exp(-beta mu) = 1/fugacity
 */
double gsl_sf_fermi_dirac(double beta, double zeta_inverse, double E);

/* Calculate the integral over the Fermi-Dirac function
   that gives the number density. Specifically,
   
   fermi_integral_1(A) = int_0^infty dy y^(1/2)/(A e^y + 1)
 */
double gsl_sf_fermi_integral_1(double A);


/* Calculate the Fermi-Dirac integral which is minus the
   logarithmic derivative of the F1 integral. Specifically

   fermi_integral_2(A) = int_0^infty dy y^(1/2)/(A e^y + 1)^2
 */
double gsl_sf_fermi_integral_2(double A);


/* Calculate the inverse fugacity given the density
   and temperature, by solving the implicit equation
   using a simple bracketing method. Specifically,
   this solves the equation

   dens_hat = g_factor/(4 Pi^2) int_0^infty dy y^(1/2)/(zeta_inv e^y + 1)

   This is normalized such if the actual density is D, then
   
   dens_hat = D  (beta/(2m))^(3/2).

   The desired precision is specified.
 */
double gsl_sf_fermi_zeta_inverse(double dens_hat, double g, double prec);


#endif /* !GSL_SF_FERMI_DIRAC_H_ */
