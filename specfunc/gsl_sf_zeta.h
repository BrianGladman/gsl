/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_ZETA_H_
#define GSL_SF_ZETA_H_


/* Riemann zeta;  zeta(x) = Sum[ k^(-x), {k,1,Infinity} ] */

int gsl_sf_zeta_e(double x, double * result);
int gsl_sf_zeta_int_e(int n, double * result);

double gsl_sf_zeta(double x);
double gsl_sf_zeta_int(int n);


/* Hurwicz zeta;  zeta(x,q) = Sum[ (k+q)^(-x), {k,0,Infinity} ] */

int gsl_sf_hzeta_e(double x, double q, double * result);

double gsl_sf_hzeta(double x, double q);


int gsl_sf_zeta_impl(double x, double * result);
int gsl_sf_hzeta_impl(double x, double q, double * result);
int gsl_sf_zeta_int_impl(int x, double * result);


#endif  /* !GSL_SF_ZETA_H_ */
