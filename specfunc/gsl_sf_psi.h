/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_PSI_H_
#define GSL_SF_PSI_H_


/* Psi(n), n > 0
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_int_impl(int n, double * result);
int     gsl_sf_psi_int_e(int n, double * result);
double  gsl_sf_psi_int(int n);


/* Psi(x), x != 0.0
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_impl(double x, double * result);
int     gsl_sf_psi_e(double x, double * result);
double  gsl_sf_psi(double x);


/* Re[Psi(1 + I y)]
 *
 * exceptions: none
 */
int     gsl_sf_psi_1piy_impl(double y, double * result);
int     gsl_sf_psi_1piy_e(double y, double * result);
double  gsl_sf_psi_1piy(double y);


#endif /* !GSL_SF_PSI_H_ */
