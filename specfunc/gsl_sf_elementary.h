/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Miscellaneous elementary functions and operations.
 */
#ifndef GSL_SF_ELEMENTARY_H_
#define GSL_SF_ELEMENTARY_H_

#include <gsl/gsl_sf_result.h>


/* Multiplication.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_multiply_impl(double x, double y, gsl_sf_result * result);
int gsl_sf_multiply_e(double x, double y, gsl_sf_result * result);


/* Multiplication of quantities with associated errors.
 */
int gsl_sf_multiply_err_impl(double x, double dx, double y, double dy, gsl_sf_result * result);
int gsl_sf_multiply_err_e(double x, double dx, double y, double dy, gsl_sf_result * result);


#endif /* !GSL_SF_ELEMENTARY_H_ */
