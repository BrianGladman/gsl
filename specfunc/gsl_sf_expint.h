/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef GSL_EXPINT_H_
#define GSL_EXPINT_H_

int gsl_sf_expint_E1_e(double x, double * result);   /* GSL_EDOM, GSL_EUNDRFLW */
int gsl_sf_expint_Ei_e(double x, double * result);   /* GSL_EDOM, GSL_EUNDRFLW */

double gsl_sf_expint_E1(double x);  /* domain, underflow */
double gsl_sf_expint_Ei(double x);  /* domain, underflow */


#endif /* !GSL_EXPINT_H_ */
