/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef GSL_SF_TRIG_H_
#define GSL_SF_TRIG_H_


/* sin(z) for complex z */
int gsl_sf_complex_sin(double zr, double zi, double * szr, double * szi);

/* cos(z) for complex z */
int gsl_sf_complex_cos(double zr, double zi, double * czr, double * czi);

int gsl_sf_complex_logsin(double zr, double zi, double * lszr, double * lszi);

/* convert polar to rectlinear coordinates */
int gsl_sf_polar_to_rect(double r, double theta, double * x, double * y);

/* convert rectilinear to polar coordinates
   return argument in range [-pi, pi]
 */
int gsl_sf_rect_to_polar(double x, double y, double * r, double * theta);

/* force an angle to lie in the range (-pi,pi] */
int gsl_sf_angle_restrict_symm(double * theta);

/* force an angle to lie in the range [0, 2pi) */
int gsl_sf_angle_restrict_pos(double * theta);


#endif /* GSL_SF_TRIG_H_ */
