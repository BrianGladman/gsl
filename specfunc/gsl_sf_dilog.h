/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_DILOG_H_
#define GSL_SF_DILOG_H_

#include <gsl_sf_result.h>


/* Real part of DiLogarithm(x), for real argument.
 * In Lewin's notation, this is Li_2(x).
 *
 *   Li_2(x) = - Re[ Integrate[ Log[1-s] / s, {s, 0, x}] ]
 *
 * Note that Im[Li_2(x)] = { 0 for x <= 1, -Pi*log(x) for x > 1 }
 */
int     gsl_sf_dilog_impl(double x, gsl_sf_result * result);
int     gsl_sf_dilog_e(double x, gsl_sf_result * result);


/* DiLogarithm(z), for complex argument z = r Exp[i theta].
 */
int gsl_sf_complex_dilog_impl(double r, double theta, gsl_sf_result * result_re, gsl_sf_result * result_im);
int gsl_sf_complex_dilog_e(double r, double theta, gsl_sf_result * result_re, gsl_sf_result * result_im);


#endif /* !GSL_SF_DILOG_H_ */
