/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_GAMMAFUNCTION_H_
#define GSL_GAMMAFUNCTION_H_


/* Log[Gamma(x)]; Lanczos method
 * Return GSL_SUCCESS on success, GSL_EDOM if x <= 0
 */
int gsl_sf_lngamma_e(double x, double * result);


/* Log[Gamma(x)]; Lanczos method */
double gsl_sf_lngamma(double);    /* domain error can occur */


/* Gamma(z) for z complex.
   Calculates:
      lnr = log|Gamma(z)|
      arg = arg(Gamma(z))  in (-Pi, Pi]
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, double * lnr, double * arg);


/* factorial(n)
 * Returns GSL_EDOM if n<=0, GSL_OVRFLW if n too large, GSL_SUCCESS on success
 */
int gsl_sf_fact_e(int n, double * result);



#endif /* !GSL_GAMMAFUNCTION_H_ */
