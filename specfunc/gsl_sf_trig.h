/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_TRIG_H_
#define GSL_SF_TRIG_H_

#include <gsl_sf_result.h>


/* sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_sin_impl(double zr, double zi, gsl_sf_result * szr, gsl_sf_result * szi);
int gsl_sf_complex_sin_e(double zr, double zi, gsl_sf_result * szr, gsl_sf_result * szi);


/* cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_cos_impl(double zr, double zi, gsl_sf_result * czr, gsl_sf_result * czi);
int gsl_sf_complex_cos_e(double zr, double zi, gsl_sf_result * czr, gsl_sf_result * czi);


/* log(sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_complex_logsin_impl(double zr, double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);
int gsl_sf_complex_logsin_e(double zr, double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);


/* sin(pi x)
 *
 * exceptions: none
 */
int gsl_sf_sin_pi_x_impl(double x, gsl_sf_result * result);
int gsl_sf_sin_pi_x_e(double x, gsl_sf_result * result);


/* log(sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnsinh_impl(double x, gsl_sf_result * result);
int gsl_sf_lnsinh_e(double x, gsl_sf_result * result);


/* log(cosh(x))
 *
 * exceptions: none
 */
int gsl_sf_lncosh_impl(double x, gsl_sf_result * result);
int gsl_sf_lncosh_e(double x, gsl_sf_result * result);


/* Convert polar to rectlinear coordinates.
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_polar_to_rect_impl(double r, double theta, gsl_sf_result * x, gsl_sf_result * y);
int gsl_sf_polar_to_rect_e(double r, double theta, gsl_sf_result * x, gsl_sf_result * y); 


/* Convert rectilinear to polar coordinates.
 * return argument in range [-pi, pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_rect_to_polar_impl(double x, double y, gsl_sf_result * r, gsl_sf_result * theta);
int gsl_sf_rect_to_polar_e(double x, double y, gsl_sf_result * r, gsl_sf_result * theta); 


/* Force an angle to lie in the range (-pi,pi].
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_symm_impl(double * theta);
int gsl_sf_angle_restrict_symm_e(double * theta);


/* Force an angle to lie in the range [0, 2pi)
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_pos_impl(double * theta);
int gsl_sf_angle_restrict_pos_e(double * theta);


/* Sin(x) for quantity with an associated error.
 */
int gsl_sf_sin_err_impl(double x, double dx, gsl_sf_result * result);
int gsl_sf_sin_err_e(double x, double dx, gsl_sf_result * result);


/* Cos(x) for quantity with an associated error.
 */
int gsl_sf_cos_err_impl(double x, double dx, gsl_sf_result * result);
int gsl_sf_cos_err_e(double x, double dx, gsl_sf_result * result);


#endif /* GSL_SF_TRIG_H_ */
