/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_trig.h"

extern int gsl_sf_complex_log_impl(double, double, double *, double *);


/*-*-*-*-*-*-*-*-*-*-*-* Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_sin_impl(double zr, double zi, double * szr, double * szi)
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
    return GSL_EOVRFLW;
  }
}

int gsl_sf_complex_logsin_impl(double zr, double zi, double * lszr, double * lszi)
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
      return GSL_EDOM;
    }
  }
  gsl_sf_angle_restrict_symm(lszi);
  return GSL_SUCCESS;
}


int gsl_sf_complex_cos_impl(double zr, double zi, double * czr, double * czi)
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
    return GSL_EOVRFLW;
  }
}


int gsl_sf_polar_to_rect_impl(double r, double theta, double * x, double * y)
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

int gsl_sf_rect_to_polar_impl(double x, double y, double * r, double * theta)
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

int gsl_sf_angle_restrict_symm_impl(double * theta, double precision)
{
  int status;
  double x;
  
  if(fabs(*theta) >  precision/GSL_MACH_EPS) {
    status = GSL_ELOSS;
  }
  else {
    status = GSL_SUCCESS;
  }
  
  x = *theta/(2.*M_PI);
  *theta = (x - floor(x)) * 2.*M_PI;
  if(*theta > M_PI) *theta -= 2.*M_PI;
  return status;
}

int gsl_sf_angle_restrict_pos_impl(double * theta, double precision)
{
  int status;
  double x;

  if(fabs(*theta) >  precision/GSL_MACH_EPS) {
    status = GSL_ELOSS;
  }
  else {
    status = GSL_SUCCESS;
  }
  
  x = *theta/(2.*M_PI);
  *theta = (x - floor(x)) * 2.*M_PI;
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_complex_sin_e(double zr, double zi, double * szr, double * szi)
{
  int status = gsl_sf_complex_sin_impl(zr, zi, szr, szi);
  
  if(status != GSL_SUCCESS) {
    char buff[128];
    sprintf(buff, "gsl_sf_complex_sin_e: n= %d", n);
    GSL_ERROR(buff, status);
  }
}

int gsl_sf_complex_logsin_e(double zr, double zi, double * lszr, double * lszi)
{
  int status = gsl_sf_complex_logsin_impl(zr, zi, lszr, lszi);
  
  if(status != GSL_SUCCESS) {
  }
}

int gsl_sf_complex_cos_e(double zr, double zi, double * czr, double * czi)
{
  int status = gsl_sf_complex_cos_impl(zr, zi, czr, czi);
  
  if(status != GSL_SUCCESS) {
  }
}


int gsl_sf_polar_to_rect_e(double r, double theta, double * x, double * y)
{
  int status = gsl_sf_polar_to_rect_impl(r, theta, x, y);
  
  if(status != GSL_SUCCESS) {
  }
}

int gsl_sf_rect_to_polar_e(double x, double y, double * r, double * theta)
{
  int status = gsl_sf_rect_to_polar_impl(x, y, r, theta);
  
  if(status != GSL_SUCCESS) {
  }
}

int gsl_sf_angle_restrict_symm_e(double * theta, double precision)
{
  int status = gsl_sf_angle_restrict_symm_impl(theta, precision);
  
  if(status != GSL_SUCCESS) {
    char buff[100];
    sprintf(buff,"gsl_sf_angle_restrict_symm_e: theta= %22.17g", *theta);
    GSL_ERROR(buff, status);
  }
}

int gsl_sf_angle_restrict_pos_e(double * theta, double precision)
{
  int status = gsl_sf_angle_restrict_pos_impl(theta, precision);
  
  if(status != GSL_SUCCESS) {
    char buff[100];
    sprintf(buff,"gsl_sf_angle_restrict_pos_e: theta= %22.17g", *theta);
    GSL_ERROR(buff, status);
  }
}
