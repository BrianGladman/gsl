/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_DAWSON_H_
#define GSL_SF_DAWSON_H_


/* Dawson's integral:
 *
 *   Exp[-x^2] Integral[ Exp[t^2], {t,0,x}]
 *
 * exceptions: GSL_EUNDRFLW;
 */
int     gsl_sf_dawson_impl(double x, gsl_sf_result * result);
int     gsl_sf_dawson_e(double x, gsl_sf_result * result);


#endif  /* !GSL_SF_DAWSON_H_ */
