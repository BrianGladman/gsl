/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_HYPERG_H_
#define GSL_SF_HYPERG_H_


int gsl_sf_hyperg_1F1_impl(double a, double b, double x, double * result);
int gsl_sf_hyperg_1F1_e(double a, double b, double x, double * result);

double gsl_sf_hyperg_1F1(double a, double b, double x);


int gsl_sf_hyperg_2F1_impl(double a, double b, double c, double x, double * result);
int gsl_sf_hyperg_2F1_conj_impl(double aR, double aI, double c, double x, double * result);


#endif  /* !GSL_SF_HYPERG_H_ */
