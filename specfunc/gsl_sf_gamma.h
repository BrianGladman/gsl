/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_GAMMAFUNCTION_H_
#define GSL_GAMMAFUNCTION_H_


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int     gsl_sf_lngamma_impl(double x, double * result);
int     gsl_sf_lngamma_e(double x, double * result);
double  gsl_sf_lngamma(double);   


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int     gsl_sf_lngamma_sgn_impl(double x, double * result_lg, double *sgn);
int     gsl_sf_lngamma_sgn_e(double x, double * result_lg, double * sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int     gsl_sf_gamma_impl(double x, double * result);
int     gsl_sf_gamma_e(double x, double * result);
double  gsl_sf_gamma(double x);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_gammastar_impl(double x, double * result);
int     gsl_sf_gammastar_e(double x, double * result);
double  gsl_sf_gammastar(double x);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int     gsl_sf_gammainv_impl(double x, double * result);
int     gsl_sf_gammainv_e(double x, double * result);
double  gsl_sf_gammainv(double x);


/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_lngamma_complex_impl(double zr, double zi, double * lnr, double * arg);
int gsl_sf_lngamma_complex_e(double zr, double zi, double * lnr, double * arg);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int     gsl_sf_fact_impl(unsigned int n, double * result);
int     gsl_sf_fact_e(unsigned int n, double * result);
double  gsl_sf_fact(unsigned int n);


/* n!! = n(n-2)(n-4) ... 
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int     gsl_sf_doublefact_impl(unsigned int n, double * result);
int     gsl_sf_doublefact_e(unsigned int n, double * result);
double  gsl_sf_doublefact(unsigned int n);


/* log(n!) 
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int     gsl_sf_lnfact_impl(unsigned int n, double * result);
int     gsl_sf_lnfact_e(unsigned int n, double * result);
double  gsl_sf_lnfact(unsigned int n);


/* log(n!!) 
 *
 * exceptions: none
 */
int     gsl_sf_lndoublefact_impl(unsigned int n, double * result);
int     gsl_sf_lndoublefact_e(unsigned int n, double * result);
double  gsl_sf_lndoublefact(unsigned int n);


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
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int     gsl_sf_lnpoch_impl(double a, double x, double * result);
int     gsl_sf_lnpoch_e(double a, double x, double * result);
double  gsl_sf_lnpoch(double a, double x);


/* Logarithm of Pochammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_impl(double a, double x, double * result, double * sgn);
int gsl_sf_lnpoch_sgn_e(double a, double x, double * result, double * sgn);


/* Pochammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int     gsl_sf_poch_impl(double a, double x, double * result);
int     gsl_sf_poch_e(double a, double x, double * result);
double  gsl_sf_poch(double a, double x);


/* Relative Pochammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int     gsl_sf_pochrel_impl(double a, double x, double * result);
int     gsl_sf_pochrel_e(double a, double x, double * result);
double  gsl_sf_pochrel(double a, double x);



/* Incomplete Gamma Function
 * gamma
 */

#endif /* !GSL_GAMMAFUNCTION_H_ */
