/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_GAMMAFUNCTION_H_
#define GSL_GAMMAFUNCTION_H_

#include <gsl_sf_result.h>


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_impl(double x, gsl_sf_result * result);
int gsl_sf_lngamma_e(double x, gsl_sf_result * result);


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_sgn_impl(double x, gsl_sf_result * result_lg, double *sgn);
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double * sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int gsl_sf_gamma_impl(double x, gsl_sf_result * result);
int gsl_sf_gamma_e(double x, gsl_sf_result * result);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gammastar_impl(double x, gsl_sf_result * result);
int gsl_sf_gammastar_e(double x, gsl_sf_result * result);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int gsl_sf_gammainv_impl(double x, gsl_sf_result * result);
int gsl_sf_gammainv_e(double x, gsl_sf_result * result);


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
int gsl_sf_lngamma_complex_impl(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


/* x^n / n!
 *
 * x >= 0.0, n >= 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_taylorcoeff_impl(int n, double x, gsl_sf_result * result);
int gsl_sf_taylorcoeff_e(int n, double x, gsl_sf_result * result);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int gsl_sf_fact_impl(unsigned int n, gsl_sf_result * result);
int gsl_sf_fact_e(unsigned int n, gsl_sf_result * result);


/* n!! = n(n-2)(n-4) ... 
 *
 * exceptions: GSL_EDOM, GSL_OVRFLW
 */
int gsl_sf_doublefact_impl(unsigned int n, gsl_sf_result * result);
int gsl_sf_doublefact_e(unsigned int n, gsl_sf_result * result);


/* log(n!) 
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int gsl_sf_lnfact_impl(unsigned int n, gsl_sf_result * result);
int gsl_sf_lnfact_e(unsigned int n, gsl_sf_result * result);


/* log(n!!) 
 *
 * exceptions: none
 */
int gsl_sf_lndoublefact_impl(unsigned int n, gsl_sf_result * result);
int gsl_sf_lndoublefact_e(unsigned int n, gsl_sf_result * result);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM 
 */
int gsl_sf_lnchoose_impl(unsigned int n, unsigned int m, gsl_sf_result * result);
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_choose_impl(unsigned int n, unsigned int m, gsl_sf_result * result);
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);


/* Logarithm of Pochammer (Apell) symbol
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_lnpoch_e(double a, double x, gsl_sf_result * result);


/* Logarithm of Pochammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_impl(double a, double x, gsl_sf_result * result, double * sgn);
int gsl_sf_lnpoch_sgn_e(double a, double x, gsl_sf_result * result, double * sgn);


/* Pochammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_poch_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_poch_e(double a, double x, gsl_sf_result * result);


/* Relative Pochammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_pochrel_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_pochrel_e(double a, double x, gsl_sf_result * result);


/* Normalized Incomplete Gamma Function
 *
 * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_Q_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_gamma_inc_Q_e(double a, double x, gsl_sf_result * result);


/* Complementary Normalized Incomplete Gamma Function
 *
 * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_P_impl(double a, double x, gsl_sf_result * result);
int gsl_sf_gamma_inc_P_e(double a, double x, gsl_sf_result * result);


/* Logarithm of Beta Function
 * Log[B(a,b)]
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnbeta_impl(double a, double b, gsl_sf_result * result);
int gsl_sf_lnbeta_e(double a, double b, gsl_sf_result * result);


/* Beta Function
 * B(a,b)
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_beta_impl(double a, double b, gsl_sf_result * result);
int gsl_sf_beta_e(double a, double b, gsl_sf_result * result);


#endif /* !GSL_GAMMAFUNCTION_H_ */
