#ifndef GSL_GAMMAFUNCTION_H_
#define GSL_GAMMAFUNCTION_H_


/* Logarithm of Gamma(x).
   Lanczos method.
 */
double gsl_sf_lngamma(double);


/* Gamma(z) for z complex.
   Calculates:
      lnr = log|Gamma(z)|
      arg = arg(Gamma(z))  in (-Pi, Pi]
 */
void gsl_sf_complex_lngamma(double zr, double zi, double * lnr, double * arg);


#endif /* GSL_GAMMAFUNCTION_H_ */
