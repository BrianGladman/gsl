/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_TRIG_H_
#define GSL_SF_TRIG_H_


/* sin(z) for complex z */
int gsl_sf_complex_sin_e(double zr, double zi, double * szr, double * szi);      /* GSL_EOVRFLW */

/* cos(z) for complex z */
int gsl_sf_complex_cos_e(double zr, double zi, double * czr, double * czi);      /* GSL_EOVRFLW */

/* log(sin(z)) */
int gsl_sf_complex_logsin_e(double zr, double zi, double * lszr, double * lszi); /* GSL_EDOM, GSL_ELOSS */

/* convert polar to rectlinear coordinates */
int gsl_sf_polar_to_rect_e(double r, double theta, double * x, double * y);    /* none */

/* convert rectilinear to polar coordinates
   return argument in range [-pi, pi]
 */
int gsl_sf_rect_to_polar_e(double x, double y, double * r, double * theta);    /* GSL_EDOM */

/* force an angle to lie in the range (-pi,pi] */
int gsl_sf_angle_restrict_symm_e(double * theta, double precision);   /* GSL_ELOSS */

/* force an angle to lie in the range [0, 2pi) */
int gsl_sf_angle_restrict_pos_e(double * theta, double precision);    /* GSL_ELOSS */


#endif /* GSL_SF_TRIG_H_ */
