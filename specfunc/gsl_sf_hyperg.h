/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_HYPERG_H_
#define GSL_SF_HYPERG_H_


/* Hypergeometric function related to Bessel functions
 * 0F1[c,x] =
 *            Gamma[c]    x^(1/2(1-c)) I_{c-1}(2 Sqrt[x])
 *            Gamma[c] (-x)^(1/2(1-c)) J_{c-1}(2 Sqrt[-x])
 *
 * exceptions:
 */
int     gsl_sf_hyperg_0F1_impl(double c, double x, double * result);
int     gsl_sf_hyperg_0F1_e(double c, double x, double * result);
double  gsl_sf_hyperg_0F1(double c, double x);


/* Confluent hypergeometric function  1F1[a,b,x] = M(a,b,x)
 *
 * exceptions:
 */
int     gsl_sf_hyperg_1F1_impl(double a, double b, double x, double * result);
int     gsl_sf_hyperg_1F1_e(double a, double b, double x, double * result);
double  gsl_sf_hyperg_1F1(double a, double b, double x);


/* Confluent hypergeometric function  U(a,b,x)
 *
 * exceptions:
 */
int     gsl_sf_hyperg_U_impl(double a, double b, double x, double * result);
int     gsl_sf_hyperg_U_e(double a, double b, double x, double * result);
double  gsl_sf_hyperg_U(double a, double b, double x);


/* Gauss hypergeometric function 2F1[a,b,c,x]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_impl(double a, double b, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_e(double a, double b, double c, double x, double * result);
double  gsl_sf_hyperg_2F1(double a, double b, double c, double x);


/* Gauss hypergeometric function
 * 2F1[aR + I aI, aR - I aI, c, x]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_conj_impl(double aR, double aI, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_conj_e(double aR, double aI, double c, double x, double * result);
double  gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x);


/* Renormalized Gauss hypergeometric function
 * 2F1[a,b,c,x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_renorm_impl(double a, double b, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_renorm_e(double a, double b, double c, double x, double * result);
double  gsl_sf_hyperg_2F1_renorm(double a, double b, double c, double x);


/* Renormalized Gauss hypergeometric function
 * 2F1[aR + I aI, aR - I aI, c, x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_conj_renorm_impl(double aR, double aI, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_conj_renorm_e(double aR, double aI, double c, double x, double * result);
double  gsl_sf_hyperg_2F1_conj_renorm(double aR, double aI, double c, double x);


/* Mysterious hypergeometric function. The series representation
 * is a divergent hypergeometric series. However, for x < 0 we
 * have 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_hyperg_2F0_series_impl(double a, double b, double x, int n_trunc, double * result);
int     gsl_sf_hyperg_2F0_impl(double a, double b, double x, double * result);
int     gsl_sf_hyperg_2F0_e(double a, double b, double x, double * result);
double  gsl_sf_hyperg_2F0(double a, double b, double x);


#endif  /* !GSL_SF_HYPERG_H_ */
