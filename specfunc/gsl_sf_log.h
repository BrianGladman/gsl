/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_LOG_H_
#define GSL_SF_LOG_H_


/* complex logarithm
 *   exp(lnr + I theta) = log(zr + I zi)
 * returns argument in [-pi,pi]
 */
void gsl_sf_complex_log(double zr, double zi, double * lnr, double * theta);


#endif /* GSL_SF_LOG_H_ */
