/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_HYPERG_H_
#define GSL_SF_HYPERG_H_


/* Hypergeometric function related to Bessel functions
 *
 * exceptions:
 */
int     gsl_sf_hyperg_0F1_impl(double c, double x, double * result);
int     gsl_sf_hyperg_0F1_e(double c, double x, double * result);
double  gsl_sf_hyperg_0F1(double c, double x);


/* Confluent hypereometric function  1F1[a,b,x] = M(a,b,x)
 *
 * exceptions:
 */
int     gsl_sf_hyperg_1F1_impl(double a, double b, double x, double * result);
int     gsl_sf_hyperg_1F1_e(double a, double b, double x, double * result);
double  gsl_sf_hyperg_1F1(double a, double b, double x);


/* Gauss hypergeometric function 2F1[a,b,c,x]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_impl(double a, double b, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_e(double a, double b, double c, double x, double * result);
double  gsl_sf_hyperg_2F1(double a, double b, double c, double x);


/* Gauss hypergeometric function 2F1[aR + I aI, aR - I aI, c, x]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_conj_impl(double aR, double aI, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_conj_e(double aR, double aI, double c, double x, double * result);
double  gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x);


/* Renormalized Gauss hypergeometric function 2F1[a,b,c,x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_impl(double a, double b, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_e(double a, double b, double c, double x, double * result);
double  gsl_sf_hyperg_2F1(double a, double b, double c, double x);


/* Gauss hypergeometric function 2F1[aR + I aI, aR - I aI, c, x] / Gamma[c]
 * |x| < 1
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F1_conj_impl(double aR, double aI, double c, double x, double * result);
int     gsl_sf_hyperg_2F1_conj_e(double aR, double aI, double c, double x, double * result);
double  gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x);


/* Mysterious hypergeometric function
 *
 * exceptions:
 */
int     gsl_sf_hyperg_2F0_impl(double a, double b, double x, double * result);
int     gsl_sf_hyperg_2F0_e(double a, double b, double x, double * result);
double  gsl_sf_hyperg_2F0(double a, double b, double x);


#endif  /* !GSL_SF_HYPERG_H_ */
