/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Miscellaneous elementary functions and operations.
 */
#ifndef GSL_SF_ELEMENTARY_H_
#define GSL_SF_ELEMENTARY_H_


/* Multiplication.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int     gsl_sf_multiply_impl(double x, double y, double * result);
int     gsl_sf_multiply_e(double x, double y, double * result);
double  gsl_sf_multiply(double x, double y);



#endif /* !GSL_SF_ELEMENTARY_H_ */
