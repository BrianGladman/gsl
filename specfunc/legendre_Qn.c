/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"


/* Evaluate f_{ell+1}/f_ell
 * f_ell := Q^{b}_{a+ell}(x)
 * x > 1
 */
static
int
legendreQ_CF1_xgt1(int ell, double a, double b, double x, double * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = ell + 1.0 + a + b;
  double b1 = (2.0*(ell+1.0+a) + 1.0) * x;
  double An = b1*Anm1 + a1*Anm2;
  double Bn = b1*Bnm1 + a1*Bnm2;
  double an, bn;
  double fn = An/Bn;

  while(n < maxiter) {
    double old_fn;
    double del;
    double lna;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    lna = ell + n + a;
    an = b*b - lna*lna;
    bn = (2.0*lna + 1.0) * x;
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;

    if(fabs(del - 1.0) < 10.0*GSL_MACH_EPS) break;
  }

  *result = fn;

  if(n == maxiter)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS; 
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_Q0_impl(const double x, double * result)
{
  if(x <= -1.0 || x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.5 * log((1.0+x)/(1.0-x));
    return GSL_SUCCESS;
  }
  else {
    *result = 0.5 * log((x+1.0)/(x-1.0));
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Q1_impl(const double x, double * result)
{
  if(x <= -1.0 || x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.5 * x * log((1.0+x)/(1.0-x)) - 1.0;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.5 * x * log((x+1.0)/(x-1.0)) - 1.0;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Q2_impl(const double x, double * result)
{
  if(x <= -1.0 || x == 1.0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(x < 1.0){
    *result = 0.25 * (3.0*x*x-1.0) * log((1.0+x)/(1.0-x)) - 1.5*x;
    return GSL_SUCCESS;
  }
  else {
    *result = 0.25 * (3.0*x*x-1.0) * log((x+1.0)/(x-1.0)) - 1.5*x;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Ql_impl(const int l, const double x, double * result)
{
  if(x <= -1.0 || x == 1.0 || l < 0) {
    *result = 0.0;
    return GSL_EDOM;
  }
  else if(l == 0) {
    return gsl_sf_legendre_Q0_impl(x, result);
  }
  else if(l == 1) {
    return gsl_sf_legendre_Q1_impl(x, result);
  }
  else if(x < 1.0){
    double Q0, Q1;
    int stat_Q0 = gsl_sf_legendre_Q0_impl(x, &Q0);
    int stat_Q1 = gsl_sf_legendre_Q1_impl(x, &Q1);
    double Qellm1 = Q0;
    double Qell   = Q1;
    double Qellp1;
    int ell;
    for(ell=1; ell<l; ell++) {
      Qellp1 = (x*(2.0*ell + 1.0) * Qell - ell * Qellm1) / (ell + 1.0);
      Qellm1 = Qell;
      Qell   = Qellp1;
    }
    *result = Qell;
    return GSL_ERROR_SELECT_2(stat_Q0, stat_Q1);
  }
  else {
    /* x > 1.0 */

    if(l > 10000) {
      /* uniform asymptotic */
      return GSL_EUNIMPL;
    }
    else {
      double rat;
      int stat_CF1  = legendreQ_CF1_xgt1(l, 0.0, 0.0, x, &rat);
      int stat_Q;
      double Qellp1 = rat * GSL_SQRT_DBL_MIN;
      double Qell   = GSL_SQRT_DBL_MIN;
      double Qellm1;
      int ell;
      for(ell=l; ell>0; ell--) {
        Qellm1 = (x * (2.0*ell + 1.0) * Qell - (ell+1.0) * Qellp1) / ell;
        Qellp1 = Qell;
        Qell   = Qellm1;
      }

      if(fabs(Qell) > fabs(Qellp1)) {
        double Q0;
        stat_Q = gsl_sf_legendre_Q0_impl(x, &Q0);
        *result = GSL_SQRT_DBL_MIN * Q0 / Qell;
      }
      else {
        double Q1;
        stat_Q = gsl_sf_legendre_Q1_impl(x, &Q1);
        *result = GSL_SQRT_DBL_MIN * Q1 / Qellp1;
      }

      return GSL_ERROR_SELECT_2(stat_Q, stat_CF1);
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_Q0_e(const double x, double * result)
{
  int status = gsl_sf_legendre_Q0_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q0_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Q1_e(const double x, double * result)
{
  int status = gsl_sf_legendre_Q1_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q1_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Q2_e(const double x, double * result)
{
  int status = gsl_sf_legendre_Q2_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Q2_e", status);
  }
  return status;
}


int
gsl_sf_legendre_Ql_e(const int l, const double x, double * result)
{
  int status = gsl_sf_legendre_Ql_impl(l, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_legendre_Ql_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_legendre_Q0(const double x)
{
  double y;
  int status = gsl_sf_legendre_Q0_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q0", status);
  }
  return y;
}


double
gsl_sf_legendre_Q1(const double x)
{
  double y;
  int status = gsl_sf_legendre_Q1_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q1", status);
  }
  return y;
}


double
gsl_sf_legendre_Q2(const double x)
{
  double y;
  int status = gsl_sf_legendre_Q2_impl(x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Q2", status);
  }
  return y;
}


double
gsl_sf_legendre_Ql(const int l, const double x)
{
  double y;
  int status = gsl_sf_legendre_Ql_impl(l, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_legendre_Ql", status);
  }
  return y;
}
