/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_PSI_H_
#define GSL_SF_PSI_H_


/* Poly-Gamma Functions
 *
 * psi(m,x) := (d/dx)^m psi(0,x) = (d/dx)^{m+1} log(gamma(x))
 */


/* Di-Gamma Function  psi(n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_int_impl(int n, double * result);
int     gsl_sf_psi_int_e(int n, double * result);
double  gsl_sf_psi_int(int n);


/* Di-Gamma Function psi(x)
 *
 * x != 0.0
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_impl(double x, double * result);
int     gsl_sf_psi_e(double x, double * result);
double  gsl_sf_psi(double x);


/* Di-Gamma Function Re[psi(1 + I y)]
 *
 * exceptions: none
 */
int     gsl_sf_psi_1piy_impl(double y, double * result);
int     gsl_sf_psi_1piy_e(double y, double * result);
double  gsl_sf_psi_1piy(double y);


/* Tri-Gamma Function psi^(1)(n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_1_int_impl(int n, double * result);
int     gsl_sf_psi_1_int_e(int n, double * result);
double  gsl_sf_psi_1_int(int n);


/* Poly-Gamma Function psi^(n)(x)
 *
 * n >= 0, x > 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_n_impl(int n, double x, double * result);
int     gsl_sf_psi_n_e(int n, double x, double * result);
double  gsl_sf_psi_n(int n, double x);


#endif /* !GSL_SF_PSI_H_ */
