/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_Q0_impl(double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.5 * log((1.0+x)/(1.0-x));
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = 0.5 * log((x+1.0)/(x-1.0));
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Q1_impl(double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.5 * x * log((1.0+x)/(1.0-x)) - 1.0;
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = 0.5 * x * log((x+1.0)/(x-1.0)) - 1.0;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Q2_impl(double x, double * result)
{
  if(x <= -1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.25 * (3.0*x*x-1.0) * log((1.0+x)/(1.0-x)) - 1.5*x;
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else {
    *result = 0.25 * (3.0*x*x-1.0) * log((x+1.0)/(x-1.0)) - 1.5*x;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_Q0_e(double x, double * result)
{
  int status = gsl_sf_legendre_Q0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q0_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Q1_e(double x, double * result)
{
  int status = gsl_sf_legendre_Q1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q1_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Q2_e(double x, double * result)
{
  int status = gsl_sf_legendre_Q2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q2_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_legendre_Q0(double x)
{
  double y;
  int status = gsl_sf_legendre_Q0_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q0", status);
  }
  return y;
}


double
gsl_sf_legendre_Q1(double x)
{
  double y;
  int status = gsl_sf_legendre_Q1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q1", status);
  }
  return y;
}


double
gsl_sf_legendre_Q2(double x)
{
  double y;
  int status = gsl_sf_legendre_Q2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q2", status);
  }
  return y;
}
