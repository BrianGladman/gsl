/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_PSI_H_
#define GSL_SF_PSI_H_


/* Psi(n)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_int_impl(int n, double * result);
int     gsl_sf_psi_int_e(int n, double * result);
double  gsl_sf_psi_int(int n);


/* Psi(x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_impl(double x, double * result);
int     gsl_sf_psi_e(double x, double * result);
double  gsl_sf_psi(double x);


#endif /* !GSL_SF_PSI_H_ */
