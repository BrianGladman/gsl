/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_DILOG_H_
#define GSL_SF_DILOG_H_


/* Real part of DiLogarithm(x), for real argument.
 * In Lewin's notation, this is Li_2(x, theta=0).
 *
 *   dilog(x) = - Re[ Integrate[ Log[1-s] / s, {s, 0, x}] ]
 */
int     gsl_sf_dilog_impl(double x, double * result);
int     gsl_sf_dilog_e(double x, double * result);
double  gsl_sf_dilog(double);


/* Real part of DiLogarithm(z), for complex argument.
 * We use the notation of Lewin again, so it is
 * expressed as a function of r=modulus(z) and c=cos(arg(z)).
 *
 *   dilogc(r, c) = -1/2 Re[ Integrate[ Log[1 - 2 c s + s^2] / s , {s,0,r}] ] 
 */
int     gsl_sf_dilogc_impl(double r, double c, double * result);
int     gsl_sf_dilogc_e(double r, double c, double * result);
double  gsl_sf_dilogc(double r, double c);


#endif /* GSL_SF_DILOG_H_ */
