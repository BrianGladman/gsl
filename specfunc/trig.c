/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_trig.h"


void gsl_sf_complex_sin(double zr, double zi, double * szr, double * szi)
{
  double ex = exp(zi);
  double sh = 0.5*(ex-1./ex);
  double ch = 0.5*(ex+1./ex);
  *szr = sin(zr)*ch;
  *szi = cos(zr)*sh;
}


void gsl_sf_complex_cos(double zr, double zi, double * czr, double * czi)
{
  double ex = exp(zi);
  double sh = 0.5*(ex-1./ex);
  double ch = 0.5*(ex+1./ex);
  *czr = cos(zr)*ch;
  *czi = -sin(zr)*sh;
}


void gsl_sf_polar_to_rect(double r, double theta, double * x, double * y)
{
  double tan_half = tan(0.5 * theta);
  double den = 1. + tan_half*tan_half;
  double cos_theta = (tan_half*tan_half - 1.) / den;
  double sin_theta = 2. * tan_half / den;
  *x = r * cos_theta;
  *y = r * sin_theta;
}


extern double hypot(double, double);

void gsl_sf_rect_to_polar(double x, double y, double * r, double * theta)
{
  *r = hypot(x, y);
  *theta = atan2(y, x);
}


void gsl_sf_angle_restrict_symm(double * theta)
{
  if(fabs(*theta) >  1.e-6 / GSL_MACH_EPS) {
    char buff[100];
    sprintf(buff,"gsl_sf_angle_restrict_symm: loss of precision for theta= %g",
	    *theta);
    GSL_ERROR(buff, GSL_EDOM);
  }
  else {
    int i_2pi = (int) floor(*theta/(2.*M_PI));
    *theta -= 2.*M_PI*i_2pi;
    if(*theta > M_PI) *theta -= 2.*M_PI;
  }
}

void gsl_sf_angle_restrict_pos(double * theta)
{
  if(fabs(*theta) > 1.e-6 / GSL_MACH_EPS) {
    char buff[100];
    sprintf(buff,"gsl_sf_angle_restrict_pos: loss of precision for theta= %g",
	    *theta);
    GSL_ERROR(buff, GSL_EDOM);
  }
  else {
    int i_2pi = (int) floor(*theta/(2.*M_PI));
    *theta -= 2.*M_PI*i_2pi;
    if(*theta > M_PI) *theta -= 2.*M_PI;
  }
}
