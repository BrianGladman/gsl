/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_bessel.h"
#include "gsl_sf_hyperg.h"

#define locEPS        (1000.0*GSL_MACH_EPS)


#if 0
static
int
hyperg_0F1_series(double c, double x, double * result, double * prec)
{
  double cn  = c;
  double n   = 1.0;
  double sum = 1.0;
  double del = 1.0;
  double abs_del = 1.0;
  double max_abs_del = 1.0;
  double err;
  int n_min = 20;
  int n_max = 200;
  int n_far = 0;
  double special_term = 0.0;
  const double rintc = floor(c + 0.5);

  /* Figure out if there is a large contribution
   * hiding far out in the sum because c is
   * near a negative integer.
   * The constant to which we compare tells how
   * many significant figures can be lost when
   * near termination if we do not do some
   * special handling.
   */
  if(c < 0.0 && fabs(c - rintc) < 0.001) {
    double ln_term;
    double ln_poch;
    double ln_fact;
    int stat_exp;
    int stat_poch;
    double sign = (x < 0.0 && GSL_IS_ODD(n_far) ? -1.0 : 1.0);
    n_far = rintc;
    gsl_sf_lnfact_impl(n_far, &ln_fact);
    stat_poch = gsl_sf_lnpoch_impl(c, n_far, &ln_poch);
    if(stat_poch != GSL_SUCCESS) {
      /* pochammer probably blew up, so there is no big contribution */
      special_term = 0.0;
      n_far = 0;
    }
    else {
      ln_term   = n_far * log(fabs(x)) - ln_fact - ln_poch;
      stat_exp  = gsl_sf_exp_impl(ln_term, &special_term);
      if(stat_exp == GSL_SUCCESS) {
        special_term *= sign;
      }
      else if(stat_exp == GSL_EUNDRFLW) {
        special_term = 0.0;
        n_far = 0;
      }
      else {
        *result = 0.0;
        return stat_exp;
      }
    }
  }

  while(abs_del/fabs(sum) > GSL_MACH_EPS || n < n_min) {
    double u, abs_u;

    if(cn == 0.0) {
      *result = 0.0;
      return GSL_EDOM;
    }
    if(n > n_max) {
      max_abs_del *= GSL_MACH_EPS;
      err     = fabs(GSL_MACH_EPS * n + max_abs_del);
      *prec   = err/(err + fabs(sum));
      *result = sum;
      if(*prec > locEPS)
        return GSL_ELOSS;
      else
        return GSL_SUCCESS;
    }
    
    u = x/(cn*n);
    abs_u = fabs(u);
    if(abs_u > 1.0 && max_abs_del > DBL_MAX/abs_u) {
      *prec   = 1.0;
      *result = sum;
      return GSL_ELOSS;
    }
    del *= u;

    sum += del;
    abs_del = fabs(del);
    max_abs_del = GSL_MAX(abs_del, max_abs_del);

    cn += 1.0;
    n  += 1.0;
  }

  /* If the sum stopped before getting to the
   * distant large contribution, then include it now.
   */
  if(n_far > n) {
    sum += special_term;
  }

  max_abs_del *= GSL_MACH_EPS;
  err     = fabs(GSL_MACH_EPS * n + max_abs_del);
  *prec   = err/(err + fabs(sum));
  *result = sum;
  if(*prec > locEPS)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}
#endif /* 0 */


/* Evaluate bessel_I(nu, x), allowing nu < 0.
 * This is fine here because we do not not allow
 * nu to be a negative integer.
 * x > 0.
 */
static
int
hyperg_0F1_bessel_I(const double nu, const double x, double * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    *result = 0.0;
    return GSL_EOVRFLW;
  }

  if(nu < 0.0) { 
    const double anu = -nu;
    const double s   = 2.0/M_PI * sin(anu*M_PI);
    const double ex  = exp(x);
    double I = 0.0;
    double K = 0.0;
    int stat_I = gsl_sf_bessel_Inu_scaled_impl(anu, x, &I);
    int stat_K = gsl_sf_bessel_Knu_scaled_impl(anu, x, &K);
    *result = ex * I + s * K / ex;
    if(stat_K != GSL_SUCCESS)
      return stat_K;
    else
      return stat_I;
  }
  else {
    const double ex  = exp(x);
    double I = 0.0;
    int stat_I = gsl_sf_bessel_Inu_scaled_impl(nu, x, &I);
    *result = ex * I;
    return stat_I;
  }
}


/* Evaluate bessel_J(nu, x), allowing nu < 0.
 * This is fine here because we do not not allow
 * nu to be a negative integer.
 * x > 0.
 */
static
int
hyperg_0F1_bessel_J(const double nu, const double x, double * result)
{
  if(nu < 0.0) { 
    const double anu = -nu;
    const double s   = sin(anu*M_PI);
    const double c   = cos(anu*M_PI);
    double J = 0.0;
    double Y = 0.0;
    int stat_J = gsl_sf_bessel_Jnu_impl(anu, x, &J);
    int stat_Y = gsl_sf_bessel_Ynu_impl(anu, x, &Y);
    *result = c * J - s * Y;
    if(stat_Y != GSL_SUCCESS)
      return stat_Y;
    else
      return stat_J;
  }
  else {
    return gsl_sf_bessel_Jnu_impl(nu, x, result);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_hyperg_0F1_impl(double c, double x, double * result)
{
  const double rintc = floor(c + 0.5);
  const int c_neg_integer = (c < 0.0 && fabs(c - rintc) < locEPS);

  if(c == 0.0 || c_neg_integer) {
    *result = 0.0;
    return GSL_EDOM;
  }

  if(x < 0.0) {
    double Jcm1;
    double lg_c, sgn;
    int stat_g = gsl_sf_lngamma_sgn_impl(c, &lg_c, &sgn);
    int stat_J = hyperg_0F1_bessel_J(c-1.0, 2.0*sqrt(-x), &Jcm1);
    if(stat_g != GSL_SUCCESS) {
      *result = 0.0;
      return stat_g;
    }
    if(Jcm1 == 0.0) {
      *result = 0.0;
      return stat_J;
    }
    else {
      double ln_pre = lg_c + log(-x)*0.5*(1.0-c);
      return gsl_sf_exp_mult_impl(ln_pre, sgn*Jcm1, result);
    }
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
    double Icm1;
    double lg_c, sgn;
    int stat_g = gsl_sf_lngamma_sgn_impl(c, &lg_c, &sgn);
    int stat_I = hyperg_0F1_bessel_I(c-1.0, 2.0*sqrt(x), &Icm1);
    if(stat_g != GSL_SUCCESS) {
      *result = 0.0;
      return stat_g;
    }
    if(Icm1 == 0.0) {
      *result = 0.0;
      return stat_I;
    }
    else {
      double ln_pre = log(x)*0.5*(1.0-c) + lg_c;
      return gsl_sf_exp_mult_impl(ln_pre, sgn*Icm1, result);
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_hyperg_0F1_e(const double c, const double x, double * result)
{
  int status = gsl_sf_hyperg_0F1_impl(c, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_0F1_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_hyperg_0F1(const double c, const double x)
{
  double y;
  int status = gsl_sf_hyperg_0F1_impl(c, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_0F1", status);
  }
  return y;
}
