/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_SYNCHROTRON_H_
#define GSL_SF_SYNCHROTRON_H_


/* First synchrotron function:
 *   synch1(x) = x Integral[ K_{5/3}(t), {t, x, Infinity}]
 */
int gsl_sf_synchrotron_1_impl(double x, double * result);


/* Second synchroton function:
 *   synch2(x) = x * K_{2/3}(x)
 */
int gsl_sf_synchrotron_2_impl(double x, double * result);


#endif  /* !GSL_SF_SYNCHROTRON_H_ */
