/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_ZETA_H_
#define GSL_SF_ZETA_H_


/* Riemann Zeta Function
 * zeta(n) = Sum[ k^(-n), {k,1,Infinity} ], n=integer
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int     gsl_sf_zeta_int_impl(int n, double * result);
int     gsl_sf_zeta_int_e(int n, double * result);
double  gsl_sf_zeta_int(int n);


/* Riemann Zeta Function
 * zeta(x) = Sum[ k^(-x), {k,1,Infinity} ]
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int     gsl_sf_zeta_impl(double x, double * result);
int     gsl_sf_zeta_e(double x, double * result);
double  gsl_sf_zeta(double x);


/* Hurwicz Zeta Function
 * zeta(x,q) = Sum[ (k+q)^(-x), {k,0,Infinity} ]
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
 */
int     gsl_sf_hzeta_impl(double x, double q, double * result);
int     gsl_sf_hzeta_e(double x, double q, double * result);
double  gsl_sf_hzeta(double x, double q);


#endif  /* !GSL_SF_ZETA_H_ */
