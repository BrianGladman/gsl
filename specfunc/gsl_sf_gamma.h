/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_GAMMAFUNCTION_H_
#define GSL_GAMMAFUNCTION_H_


/* Log[Gamma(x)]; Lanczos method
 * Return GSL_SUCCESS on success, GSL_EDOM if x <= 0
 */
int gsl_sf_lngamma_e(double x, double * result);  /* GSL_EDOM */

double gsl_sf_lngamma(double);    /* domain error can occur */


/* Gamma(z) for z complex.
   Calculates:
      lnr = log|Gamma(z)|
      arg = arg(Gamma(z))  in (-Pi, Pi]
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, double * lnr, double * arg);


/* n!
 * Returns GSL_EDOM if n<=0, GSL_OVRFLW if n too large, GSL_SUCCESS on success
 */
int gsl_sf_fact_e(int n, double * result);

/* log(n!) 
 * Faster then ln(Gamma(n+1)) for n < 170; defers for larger n.
 */
int gsl_sf_lnfact_e(int n, double * result);  /* GSL_EDOM */

double gsl_sf_lnfact(int n);  /* domain */


#endif /* !GSL_GAMMAFUNCTION_H_ */
