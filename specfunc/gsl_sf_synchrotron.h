/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_SYNCHROTRON_H_
#define GSL_SF_SYNCHROTRON_H_


/* First synchrotron function:
 *   synchrotron_1(x) = x Integral[ K_{5/3}(t), {t, x, Infinity}]
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_synchrotron_1_impl(double x, double * result);
int     gsl_sf_synchrotron_1_e(double x, double * result);
double  gsl_sf_synchrotron_1(double x);


/* Second synchroton function:
 *   synchrotron_2(x) = x * K_{2/3}(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_synchrotron_2_impl(double x, double * result);
int     gsl_sf_synchrotron_2_e(double x, double * result);
double  gsl_sf_synchrotron_2(double x);


#endif  /* !GSL_SF_SYNCHROTRON_H_ */
