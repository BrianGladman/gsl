#ifndef GSL_DILOG_H_
#define GSL_DILOG_H_


/* Real part of DiLogarithm(x), for real argument.
 * In Lewin's notation, this is Li_2(x, theta=0).
 *
 *   dilog(x) = - Re[ Integrate[ Log[1-s] / s, {s, 0, x}] ]
 */
double gsl_sf_dilog(double);


/* Real part of DiLogarithm(z), for complex argument.
 * We use the notation of Lewin again, so it is
 * expressed as a function of the modulus and cos(phase).
 *
 *   dilog_c(r, cos) = -1/2 Re[ Integrate[ Log[1 - 2 cos s + s^2] / s , {s,0,r}] ] 
 */
double gsl_sf_dilog_c(double r, double cos_theta);


#endif /* GSL_DILOG_H_ */
