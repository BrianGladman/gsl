/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef GSL_EXPINT_H_
#define GSL_EXPINT_H_

int gsl_sf_expint_E1_e(double x, double * result);   /* GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW */
int gsl_sf_expint_Ei_e(double x, double * result);   /* GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW */

double gsl_sf_expint_E1(double x);  /* domain, overflow, underflow */
double gsl_sf_expint_Ei(double x);  /* domain, overflow, underflow */


int gsl_sf_Shi_e(double x, double * result);  /* GSL_EOVRFLW, GSL_EUNDRFLW */
int gsl_sf_Chi_e(double x, double * result);  /* GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW */

double gsl_sf_Shi(double x);  /* overflow, underflow */
double gsl_sf_Chi(double x);  /* domain, overflow, underflow */


int gsl_sf_Si_e(double x, double * result);
int gsl_sf_Ci_e(double x, double * result);  /* GSL_EDOM */

double gsl_sf_Si(double x);
double gsl_sf_Ci(double x);  /* domain */


#endif /* !GSL_EXPINT_H_ */
