/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_TRIG_H_
#define GSL_SF_TRIG_H_


/* sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_sin_impl(double zr, double zi, double * szr, double * szi);
int gsl_sf_complex_sin_e(double zr, double zi, double * szr, double * szi);


/* cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_cos_impl(double zr, double zi, double * czr, double * czi);
int gsl_sf_complex_cos_e(double zr, double zi, double * czr, double * czi);


/* log(sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_complex_logsin_impl(double zr, double zi, double * lszr, double * lszi);
int gsl_sf_complex_logsin_e(double zr, double zi, double * lszr, double * lszi);


/* log(sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int     gsl_sf_lnsinh_impl(double x, double * result);
int     gsl_sf_lnsinh_e(double x, double * result);
double  gsl_sf_lnsinh(double x);


/* log(cosh(x))
 *
 * exceptions: none
 */
int     gsl_sf_lncosh_impl(double x, double * result);
int     gsl_sf_lncosh_e(double x, double * result);
double  gsl_sf_lncosh(double x);


/* Convert polar to rectlinear coordinates.
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_polar_to_rect_impl(double r, double theta, double * x, double * y);
int gsl_sf_polar_to_rect_e(double r, double theta, double * x, double * y); 


/* Convert rectilinear to polar coordinates.
 * return argument in range [-pi, pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_rect_to_polar_impl(double x, double y, double * r, double * theta);
int gsl_sf_rect_to_polar_e(double x, double y, double * r, double * theta); 


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


double gsl_sf_sin_pi_x(double x);


#endif /* GSL_SF_TRIG_H_ */
