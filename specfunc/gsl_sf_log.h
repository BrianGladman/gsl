/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_LOG_H_
#define GSL_SF_LOG_H_


/* complex logarithm
 *   exp(lnr + I theta) = zr + I zi
 * returns argument in [-pi,pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_impl(double zr, double zi, double * lnr, double * theta);
int gsl_sf_complex_log_e(double zr, double zi, double * lnr, double * theta);


/* log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_log_1plusx_impl(double x, double * result);
int     gsl_sf_log_1plusx_e(double x, double * result);
double  gsl_sf_log_1plusx(double x);


#endif /* GSL_SF_LOG_H_ */
