/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_TRIG_H
#define GSL_SF_TRIG_H

#include <gsl_sf_result.h>


/* Sin(x) with GSL semantics. This is actually important
 * because we want to control the error estimate, and trying
 * to guess the error for the standard library implementation
 * every time it is used would be a little goofy.
 */
int gsl_sf_sin_impl(double x, gsl_sf_result * result);
int gsl_sf_sin_e(double x, gsl_sf_result * result);


/* Cos(x) with GSL semantics.
 */
int gsl_sf_cos_impl(double x, gsl_sf_result * result);
int gsl_sf_cos_e(double x, gsl_sf_result * result);


/* Hypot(x,y) with GSL semantics.
 */
int gsl_sf_hypot_impl(double x, double y, gsl_sf_result * result);
int gsl_sf_hypot_e(double x, double y, gsl_sf_result * result);


/* Sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_sin_impl(double zr, double zi, gsl_sf_result * szr, gsl_sf_result * szi);
int gsl_sf_complex_sin_e(double zr, double zi, gsl_sf_result * szr, gsl_sf_result * szi);


/* Cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_cos_impl(double zr, double zi, gsl_sf_result * czr, gsl_sf_result * czi);
int gsl_sf_complex_cos_e(double zr, double zi, gsl_sf_result * czr, gsl_sf_result * czi);


/* Log(Sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_complex_logsin_impl(double zr, double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);
int gsl_sf_complex_logsin_e(double zr, double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);


/* Sinc(x) = sin(pi x) / (pi x)
 *
 * exceptions: none
 */
int gsl_sf_sinc_impl(double x, gsl_sf_result * result);
int gsl_sf_sinc_e(double x, gsl_sf_result * result);


/* Log(Sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnsinh_impl(double x, gsl_sf_result * result);
int gsl_sf_lnsinh_e(double x, gsl_sf_result * result);


/* Log(Cosh(x))
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


/* Sin(x) for quantity with an associated error.
 */
int gsl_sf_sin_err_impl(double x, double dx, gsl_sf_result * result);
int gsl_sf_sin_err_e(double x, double dx, gsl_sf_result * result);


/* Cos(x) for quantity with an associated error.
 */
int gsl_sf_cos_err_impl(double x, double dx, gsl_sf_result * result);
int gsl_sf_cos_err_e(double x, double dx, gsl_sf_result * result);


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


int gsl_sf_angle_restrict_symm_err_impl(double theta, gsl_sf_result * result);

int gsl_sf_angle_restrict_pos_err_impl(double theta, gsl_sf_result * result);


#endif /* !GSL_SF_TRIG_H */
