/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_TRANSPORT_H_
#define GSL_SF_TRANSPORT_H_


/* Transport function:
 *   J(n,x) := Integral[ t^n e^t /(e^t - 1)^2, {t,0,x}]
 */

/* J(2,x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_transport_2_impl(double x, double * result);
int     gsl_sf_transport_2_e(double x, double * result);
double  gsl_sf_transport_2(double x);


/* J(3,x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_transport_3_impl(double x, double * result);
int     gsl_sf_transport_3_e(double x, double * result);
double  gsl_sf_transport_3(double x);


#endif  /* !GSL_SF_TRANSPORT_H_ */
