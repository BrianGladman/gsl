/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_trig.h"

extern int gsl_sf_complex_log(double, double, double *, double *);


int gsl_sf_complex_sin(double zr, double zi, double * szr, double * szi)
{
  if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double sh = 0.5*(ex-1./ex);
    double ch = 0.5*(ex+1./ex);
    *szr = sin(zr)*ch;
    *szi = cos(zr)*sh;
    return GSL_SUCCESS;
  }
  else {
    GSL_ERROR("gsl_sf_complex_sin: Im(z) too large", GSL_EOVRFLW);
  }
}

int gsl_sf_complex_logsin(double zr, double zi, double * lszr, double * lszi)
{
  if(zi > 60.) {
    *lszr = -M_LN2 + zi;
    *lszi =  M_PI*0.5 - zr;
  }
  else if(zi < -60.) {
    *lszr = -M_LN2 - zi;
    *lszi = -0.5*M_PI + zr; 
  }
  else {
    double sin_r, sin_i;
    int status;
    gsl_sf_complex_sin(zr, zi, &sin_r, &sin_i);
    status = gsl_sf_complex_log(sin_r, sin_i, lszr, lszi);
    if(status == GSL_EDOM) {
      GSL_ERROR("gsl_sf_complex_logsin: sin(z) = 0.0", GSL_EDOM);
    }
  }
  gsl_sf_angle_restrict_symm(lszi);
  return GSL_SUCCESS;
}


int gsl_sf_complex_cos(double zr, double zi, double * czr, double * czi)
{
  if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double sh = 0.5*(ex-1./ex);
    double ch = 0.5*(ex+1./ex);
    *czr =  cos(zr)*ch;
    *czi = -sin(zr)*sh;
    return GSL_SUCCESS;
  }
  else {
    GSL_ERROR("gsl_sf_complex_cos: Im(z) too large", GSL_EOVRFLW);
  }
}


int gsl_sf_polar_to_rect(double r, double theta, double * x, double * y)
{
  if(fabs(theta) == M_PI) {
    *x = r;
    *y = 0.;
  }
  else {
    double tan_half = tan(0.5 * theta);
    double den = 1. + tan_half*tan_half;
    double cos_theta = (tan_half*tan_half - 1.) / den;
    double sin_theta = 2. * tan_half / den;
    *x = r * cos_theta;
    *y = r * sin_theta;
  }
  return GSL_SUCCESS;
}


extern double hypot(double, double);

int gsl_sf_rect_to_polar(double x, double y, double * r, double * theta)
{
  *r = hypot(x, y);
  if(*r > 0.) {
    *theta = atan2(y, x);
  }
  else {
    *theta = 0.;
  }
  return GSL_SUCCESS;
}


int gsl_sf_angle_restrict_symm(double * theta)
{
  if(fabs(*theta) >  1.e-6 / GSL_MACH_EPS) {
    char buff[100];
    sprintf(buff,"gsl_sf_angle_restrict_symm: loss of precision for theta= %g",
	    *theta);
    GSL_ERROR(buff, GSL_ELOSS);
  }
  else {
    int i_2pi = (int) floor(*theta/(2.*M_PI));
    *theta -= 2.*M_PI*i_2pi;
    if(*theta > M_PI) *theta -= 2.*M_PI;
    return GSL_SUCCESS;
  }
}

int gsl_sf_angle_restrict_pos(double * theta)
{
  if(fabs(*theta) > 1.e-6 / GSL_MACH_EPS) {
    char buff[100];
    sprintf(buff,"gsl_sf_angle_restrict_pos: loss of precision for theta= %g",
	    *theta);
    GSL_ERROR(buff, GSL_ELOSS);
  }
  else {
    int i_2pi = (int) floor(*theta/(2.*M_PI));
    *theta -= 2.*M_PI*i_2pi;
    if(*theta > M_PI) *theta -= 2.*M_PI;
    return GSL_SUCCESS;
  }
}
