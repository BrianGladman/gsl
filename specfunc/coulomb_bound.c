/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_laguerre.h"
#include "gsl_sf_coulomb.h"


/* normalization for hydrogenic wave functions */
static
double
R_norm(const int n, const int l, const double Z)
{
  double A = 2.0*Z/n;
  double term1 = A*A*A /(2.0*n);
  double ln_a, ln_b;
  gsl_sf_lnfact_impl(n+l, &ln_a);
  gsl_sf_lnfact_impl(n-l-1, &ln_b);
  return sqrt(term1) * exp(-0.5*(ln_a-ln_b));
}


int
gsl_sf_hydrogenicR_1_impl(const double Z, const double r, double * result)
{
  if(Z > 0.0 && r >= 0.0) {
    double A = 2.0*Z;
    double norm = A*sqrt(Z);
    double ea = exp(-Z*r);
    *result = norm*ea;
    return ( *result > 0.0 ? GSL_SUCCESS : GSL_EUNDRFLW );
  }
  else {
    return GSL_EDOM;
  }
}

int
gsl_sf_hydrogenicR_2_impl(const int l, const double Z, const double r, double * result)
{
  if(Z > 0.0 && r >= 0.0) {
    double term1 = 0.25*Z*Z*Z;
    double ea = exp(-Z*r);
    if(l == 0) {
      double pp = 2.0-Z*r;
      double norm = sqrt(term1/2.0);
      *result = norm * ea * pp;
      return GSL_SUCCESS;
    }
    else if(l == 1) {
      double pp = Z*r;
      double norm = sqrt(term1/6.0);
      *result = norm * ea * pp;
      return GSL_SUCCESS;
    }
    else {
      return GSL_EDOM;
    }
  }
  else {
    return GSL_EDOM;
  }
}

int
gsl_sf_hydrogenicR_impl(const int n, const int l, const double Z, const double r, double * result)
{
  if(n < 1 || l > n-1 || Z <= 0.0) {
    return GSL_EDOM;
  }
  else {
    double A = 2.0*Z/n;
    double norm = R_norm(n, l, Z);
    double rho = A*r;
    double ea = exp(-0.5*rho);
    double pp = gsl_sf_pow_int(rho, l);
    double lag;
    gsl_sf_laguerre_n_impl(n-l-1, 2*l+1, rho, &lag);
    *result = norm * ea * pp * lag;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_hydrogenicR_1_e(const double Z, const double r, double * result)
{
  int status = gsl_sf_hydrogenicR_1_impl(Z, r, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hydrogenicR_1_e", status);
  }
  return status;
}

int
gsl_sf_hydrogenicR_2_e(const int l, const double Z, const double r, double * result)
{
  int status = gsl_sf_hydrogenicR_2_impl(l, Z, r, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hydrogenicR_2_e", status);
  }
  return status;
}

int
gsl_sf_hydrogenicR_e(const int n, const int l, const double Z, const double r, double * result)
{
  int status = gsl_sf_hydrogenicR_impl(n, l, Z, r, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hydrogenicR_e", status);
  }
  return status;
}


double
gsl_sf_hydrogenicR_1(const double Z, const double r)
{
  double y;
  int status = gsl_sf_hydrogenicR_1_impl(Z, r, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hydrogenicR_1", status);
  }
  return y;
}

double
gsl_sf_hydrogenicR_2(const int l, const double Z, const double r)
{
  double y;
  int status = gsl_sf_hydrogenicR_2_impl(l, Z, r, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hydrogenicR_2", status);
  }
  return y;
}

double
gsl_sf_hydrogenicR(const int n, const int l, const double Z, const double r)
{
  double y;
  int status = gsl_sf_hydrogenicR_impl(n, l, Z, r, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hydrogenicR", status);
  }
  return y;
}
