/* Author:  G. Jungman
 * RCS:     $Id$
 */
#ifndef TRIG_IMPL_H_
#define TRIG_IMPL_H_

int gsl_sf_complex_sin_impl(double zr, double zi, double * szr, double * szi);
int gsl_sf_complex_logsin_impl(double zr, double zi, double * lszr, double * lszi);

int gsl_sf_complex_cos_impl(double zr, double zi, double * czr, double * czi);

int gsl_sf_polar_to_rect_impl(double r, double theta, double * x, double * y);
int gsl_sf_rect_to_polar_impl(double x, double y, double * r, double * theta);

int gsl_sf_angle_restrict_symm_impl(double * theta, double precision);
int gsl_sf_angle_restrict_pos_impl(double * theta, double precision);


#endif /* !TRIG_IMPL_H_ */
