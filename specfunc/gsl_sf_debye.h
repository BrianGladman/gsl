/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef __GSL_SF_DEBYE_H__
#define __GSL_SF_DEBYE_H__

#include <gsl/gsl_sf_result.h>


/* D_n(x) := n/x^n Integrate[t^n/(e^t - 1), {t,0,x}] */

/* D_1(x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_debye_1_impl(double x, gsl_sf_result * result);
int     gsl_sf_debye_1_e(double x, gsl_sf_result * result);


/* D_2(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_2_impl(double x, gsl_sf_result * result);
int     gsl_sf_debye_2_e(double x, gsl_sf_result * result);


/* D_3(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_3_impl(double x, gsl_sf_result * result);
int     gsl_sf_debye_3_e(double x, gsl_sf_result * result);


/* D_4(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_4_impl(double x, gsl_sf_result * result);
int     gsl_sf_debye_4_e(double x, gsl_sf_result * result);


#endif /* __GSL_SF_DEBYE_H__ */
