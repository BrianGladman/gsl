/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"


int gsl_sf_exp_impl(double x, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }
  else if(x < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = exp(x);
    return GSL_SUCCESS;
  }
}

int gsl_sf_expm1_impl(double x, double * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    *result = -1.0;
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    *result = exp(x) - 1.0;
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    *result = x * (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
    return GSL_SUCCESS;
  } 
  else if(x < GSL_LOG_DBL_MAX) {
    *result = exp(x) - 1.0;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
}


int gsl_sf_exprel_impl(double x, double * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    *result = -1.0/x;
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    *result = (exp(x) - 1.0)/x;
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    *result = (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
    return GSL_SUCCESS;
  } 
  else if(x < GSL_LOG_DBL_MAX) {
    *result = (exp(x) - 1.0)/x;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
}


int gsl_sf_exp_e(double x, double * result)
{
  int status = gsl_sf_exp_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exp_e", status);
  }
  return status;
}

int gsl_sf_expm1_e(double x, double * result)
{
  int status = gsl_sf_expm1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_expm1_e", status);
  }
  return status;
}


int gsl_sf_exprel_e(double x, double * result)
{
  int status = gsl_sf_exprel_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_exprel_e", status);
  }
  return status;
}


double gsl_sf_exp(double x)
{
  double y;
  int status = gsl_sf_exp_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exp", status);
  }
  return y;
}


double gsl_sf_expm1(double x)
{
  double y;
  int status = gsl_sf_expm1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_expm1", status);
  }
  return y;
}


double gsl_sf_exprel(double x)
{
  double y;
  int status = gsl_sf_exprel_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("  gsl_sf_exprel", status);
  }
  return y;
}
