/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_TRANSPORT_H_
#define GSL_SF_TRANSPORT_H_

#include <gsl/gsl_sf_result.h>


/* Transport function:
 *   J(n,x) := Integral[ t^n e^t /(e^t - 1)^2, {t,0,x}]
 */

/* J(2,x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_transport_2_impl(double x, gsl_sf_result * result);
int     gsl_sf_transport_2_e(double x, gsl_sf_result * result);


/* J(3,x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_transport_3_impl(double x, gsl_sf_result * result);
int     gsl_sf_transport_3_e(double x, gsl_sf_result * result);


/* J(4,x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_transport_4_impl(double x, gsl_sf_result * result);
int     gsl_sf_transport_4_e(double x, gsl_sf_result * result);


/* J(5,x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_transport_5_impl(double x, gsl_sf_result * result);
int     gsl_sf_transport_5_e(double x, gsl_sf_result * result);


#endif  /* !GSL_SF_TRANSPORT_H_ */
