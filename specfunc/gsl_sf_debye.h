/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_DEBYE_H_
#define GSL_SF_DEBYE_H_


/* D_n(x) := n/x^n Integrate[t^n/(e^t - 1), {t,0,x}] */

/* D_1(x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_debye_1_impl(double x, double * result);
int     gsl_sf_debye_1_e(double x, double * result);
double  gsl_sf_debye_1(double x);


/* D_2(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_2_impl(double x, double * result);
int     gsl_sf_debye_2_e(double x, double * result);
double  gsl_sf_debye_2(double x);


/* D_3(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_3_impl(double x, double * result);
int     gsl_sf_debye_3_e(double x, double * result);
double  gsl_sf_debye_3(double x);


/* D_4(x)
 *
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int     gsl_sf_debye_4_impl(double x, double * result);
int     gsl_sf_debye_4_e(double x, double * result);
double  gsl_sf_debye_4(double x);


#endif /* !GSL_SF_DEBYE_H_ */
