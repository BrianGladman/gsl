/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_legendre.h"


int
gsl_sf_legendre_H3d_0_impl(const double lambda, const double eta, double * result)
{
  if(eta < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(eta == 0.0 || lambda == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
    if(eta + log(fabs(lambda)) > -GSL_LOG_DBL_MIN) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else if(eta > -0.5*GSL_LOG_MACH_EPS) {
      *result = 2.0 * sin(lambda*eta)/lambda * exp(-eta);
      return GSL_SUCCESS;
    }
    else {
      *result = sin(lambda*eta)/(lambda*sinh(eta));
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_H3d_1_impl(const double lambda, const double eta, double * result)
{
  const double xi    = fabs(eta*lambda);
  const double lsqp1 = lambda*lambda + 1.0;

  if(eta < 0.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(eta == 0.0 || lambda == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else if(eta < GSL_ROOT4_MACH_EPS && xi < GSL_ROOT4_MACH_EPS) {
    *result = eta * sqrt(lsqp1)/3.0 * (1.0 - (2.0 + 3.0*lambda*lambda)*eta*eta/30.0);
    if(*result == 0.0)
      return GSL_EUNDRFLW;
    else
      return GSL_SUCCESS;
  }
  else {
    double sin_term;
    double cos_term;
    double coth_term;
    if(xi < GSL_ROOT4_MACH_EPS) {
      sin_term = 1.0 - xi*xi/6.0;
      cos_term = 1.0 - 0.5*xi*xi;
    }
    else {
      sin_term = sin(xi)/xi;
      cos_term = cos(xi);
    }
    if(eta < GSL_ROOT4_MACH_EPS) {
      coth_term = 1.0 + eta*eta/3.0 * (1.0 - eta*eta/15.0);
    }
    else {
      coth_term = eta/tanh(eta);
    }
    *result = 1.0/sqrt(lsqp1)/eta * (coth_term * sin_term - cos_term);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_hyper_array_impl(int lmax, double lambda, double x, double * result, double * harvest)
{
  double X = 1./tanh(x);
  double y2, y1, y0;
  int ell;

/*
  gsl_sf_hyper_0_impl(lambda, x, &y2);
  gsl_sf_hyper_1_impl(lambda, x, &y1);
*/
  harvest[0] = y2;
  harvest[1] = y1;

  for(ell=2; ell<=lmax; ell++) {
    double a = sqrt(lambda*lambda + ell*ell);
    double b = sqrt(lambda*lambda + (ell-1)*(ell-1));
    y0 = ((2*ell-1)*X*y1 - b*y2) / a;
    y2 = y1;
    y1 = y0;
    harvest[ell] = y0;
  }

  *result = y0;
}



int
gsl_sf_legendre_H3d_0_e(const double lambda, const double eta, double * result)
{
  int status = gsl_sf_legendre_H3d_0_impl(lambda, eta, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_H3d_0_e", status);
  }
  return status;
}


int
gsl_sf_legendre_H3d_1_e(const double lambda, const double eta, double * result)
{
  int status = gsl_sf_legendre_H3d_1_impl(lambda, eta, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_H3d_1_e", status);
  }
  return status;
}


double
gsl_sf_legendre_H3d_0(const double lambda, const double eta)
{
  double y;
  int status = gsl_sf_legendre_H3d_0_impl(lambda, eta, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_H3d_0", status);
  }
  return y;
}


double
gsl_sf_legendre_H3d_1(const double lambda, const double eta)
{
  double y;
  int status = gsl_sf_legendre_H3d_1_impl(lambda, eta, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_H3d_1", status);
  }
  return y;
}
