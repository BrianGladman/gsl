/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_HYPERG_H_
#define GSL_SF_HYPERG_H_


/* Confluent hypereometric function */

int gsl_sf_hyperg_1F1_impl(double a, double b, double x, double * result);
int gsl_sf_hyperg_1F1_e(double a, double b, double x, double * result);

double gsl_sf_hyperg_1F1(double a, double b, double x);


/* Gauss hypergeometric function */

int gsl_sf_hyperg_2F1_impl(double a, double b, double c, double x, double * result);
int gsl_sf_hyperg_2F1_conj_impl(double aR, double aI, double c, double x, double * result);

int gsl_sf_hyperg_2F1_e(double a, double b, double c, double x, double * result);
int gsl_sf_hyperg_2F1_conj_e(double aR, double aI, double c, double x, double * result);

double gsl_sf_hyperg_2F1(double a, double b, double c, double x);
double gsl_sf_hyperg_2F1_conj(double aR, double aI, double c, double x);


/* Mysterious hypergeometric function */

int gsl_sf_hyperg_2F0_impl(double a, double b, double x, double * result);
int gsl_sf_hyperg_2F0_e(double a, double b, double x, double * result);

double gsl_sf_hyperg_2F0(double a, double b, double x);


#endif  /* !GSL_SF_HYPERG_H_ */
