/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_GAMMAFUNCTION_H_
#define GSL_GAMMAFUNCTION_H_


/* Log[Gamma(x)]; Lanczos method
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_lngamma_impl(double x, double * result);
int     gsl_sf_lngamma_e(double x, double * result);
double  gsl_sf_lngamma(double);   


/* Gamma(z) for z complex.
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lngamma_complex_impl(double zr, double zi, double * lnr, double * arg);
int gsl_sf_lngamma_complex_e(double zr, double zi, double * lnr, double * arg);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int gsl_sf_fact_impl(unsigned int n, double * result);
int gsl_sf_fact_e(unsigned int n, double * result);


/* n!! = n(n-2)(n-4) ... 
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int gsl_sf_doublefact_impl(unsigned int n, double * result);
int gsl_sf_doublefact_e(unsigned int n, double * result);


/* log(n!) 
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_lnfact_impl(unsigned int n, double * result);
int     gsl_sf_lnfact_e(unsigned int n, double * result);
double  gsl_sf_lnfact(unsigned int n);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM 
 */
int     gsl_sf_lnchoose_impl(unsigned int n, unsigned int m, double * result);
int     gsl_sf_lnchoose_e(unsigned int n, unsigned int m, double * result);
double  gsl_sf_lnchoose(unsigned int n, unsigned int m);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int     gsl_sf_choose_impl(unsigned int n, unsigned int m, double * result);
int     gsl_sf_choose_e(unsigned int n, unsigned int m, double * result);
double  gsl_sf_choose(unsigned int n, unsigned int m);


/* Logarithm of Pochammer (Apell) symbol
 *   log( (a)_n )
 *   where (a)_n := Gamma[a + n]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int     gsl_sf_lnpoch_impl(double a, int n, double * result);
int     gsl_sf_lnpoch_e(double a, int n, double * result);
double  gsl_sf_lnpoch(double a, int n);
 
 
#endif /* !GSL_GAMMAFUNCTION_H_ */
