/* Author: G. Jungman
 * RCS: $Id$
 */
#ifndef GSL_ELLINT_H_
#define GSL_ELLINT_H_


int gsl_sf_ellint_RC_impl(double x, double y, double errtol, double * result);
int gsl_sf_ellint_RD_impl(double x, double y, double z, double errtol, double * result);
int gsl_sf_ellint_RF_impl(double x, double y, double z, double errtol, double * result);
int gsl_sf_ellint_RJ_impl(double x, double y, double z, double p, double errtol, double * result);


#endif  /* !GSL_ELLINT_H_ */
