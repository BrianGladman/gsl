/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_laguerre.h"


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* based on the large 2b-4a asymptotic for 1F1
 * [Abramowitz+Stegun, 13.5.21]
 */
static
int
laguerre_large_n(const int n, const double alpha, const double x,
                 gsl_sf_result * result)
{
  const double a = -n;
  const double b = alpha + 1.0;
  const double eta    = 2.0*b - 4.0*a;
  const double cos2th = x/eta;
  const double sin2th = 1.0 - cos2th;
  const double th = acos(sqrt(cos2th));
  const double pre_h  = 0.25*M_PI*M_PI*eta*eta*cos2th*sin2th;
  gsl_sf_result lg_b;
  gsl_sf_result lnfact;
  int stat_lg = gsl_sf_lngamma_impl(b+n, &lg_b);
  int stat_lf = gsl_sf_lnfact_impl(n, &lnfact);
  double pre_term1 = 0.5*(1.0-b)*log(0.25*x*eta);
  double pre_term2 = 0.25*log(pre_h);
  double lnpre_val = lg_b.val - lnfact.val + 0.5*x + pre_term1 - pre_term2;
  double lnpre_err = lg_b.err + lnfact.err + GSL_DBL_EPSILON * (fabs(pre_term1)+fabs(pre_term2));
  double ser_term1 = sin(a*M_PI);
  double ser_term2 = sin(0.25*eta*(2.0*th - sin(2.0*th)) + 0.25*M_PI);
  double ser_val = ser_term1 + ser_term2;
  double ser_err = GSL_DBL_EPSILON * (fabs(ser_term1) + fabs(ser_term2));
  int stat_e = gsl_sf_exp_mult_err_impl(lnpre_val, lnpre_err, ser_val, ser_err, result);
  result->err += 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
  return GSL_ERROR_SELECT_3(stat_e, stat_lf, stat_lg);
}


/* Evaluate polynomial based on confluent hypergeometric representation.
 *
 * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
 *
 * assumes n > 0
 */
static
int
laguerre_n_cp(const int n, const double a, const double x, gsl_sf_result * result)
{
  gsl_sf_result lnfact;
  gsl_sf_result lnpoch;
  int stat_f = gsl_sf_lnfact_impl(n, &lnfact);
  int stat_p = gsl_sf_lnpoch_impl(a+1.0, n, &lnpoch);
  double poly_1F1_val = 1.0;
  double poly_1F1_err = 0.0;
  int stat_e;
  int k;

  double lnpre_val = lnpoch.val - lnfact.val;
  double lnpre_err = lnpoch.err + lnfact.err + GSL_DBL_EPSILON * fabs(lnpre_val);

  for(k=n-1; k>=0; k--) {
    double t = (-n+k)/(a+1.0+k) * (x/(k+1));
    double r = t + 1.0/poly_1F1_val;
    if(r > 0.9*GSL_DBL_MAX/poly_1F1_val) {
      result->val = 0.0; /* FIXME: should be Inf */
      result->err = 0.0;
      return GSL_EOVRFLW;
    }
    else {
      /* Collect the Horner terms. */
      poly_1F1_val  = 1.0 + t * poly_1F1_val;
      poly_1F1_err += GSL_DBL_EPSILON + fabs(t) * poly_1F1_err;
    }
  }

  stat_e = gsl_sf_exp_mult_err_impl(lnpre_val, lnpre_err,
                                    poly_1F1_val, poly_1F1_err,
                                    result);

  return GSL_ERROR_SELECT_3(stat_e, stat_f, stat_p);
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_laguerre_1_impl(const double a, const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    result->val = 1.0 + a - x;
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(a) + fabs(x));
    return GSL_SUCCESS;
  }
}

int
gsl_sf_laguerre_2_impl(const double a, const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    double c0 = 0.5 * (2.0+a)*(1.0+a);
    double c1 = -(2.0+a);
    double c2 = -0.5/(2.0+a);
    result->val  = c0 + c1*x*(1.0 + c2*x);
    result->err  = GSL_DBL_EPSILON * (fabs(c0) + 2.0 * fabs(c1*x) * (1.0 + 2.0 * fabs(c2*x)));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_laguerre_3_impl(const double a, const double x, gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else {
    double c0 = (3.0+a)*(2.0+a)*(1.0+a) / 6.0;
    double c1 = -c0 * 3.0 / (1.0+a);
    double c2 = -1.0/(2.0+a);
    double c3 = -1.0/(3.0*(3.0+a));
    result->val  = c0 + c1*x*(1.0 + c2*x*(1.0 + c3*x));
    result->err  = 1.0 + 2.0 * fabs(c3*x);
    result->err  = 1.0 + 2.0 * fabs(c2*x) * result->err;
    result->err  = GSL_DBL_EPSILON * (fabs(c0) + 2.0 * fabs(c1*x) * result->err);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_laguerre_n_impl(const int n, const double a, const double x,
                           gsl_sf_result * result)
{
  if(result == 0) {
    return GSL_EFAULT;
  }
  else if(n < 0 || a <= -1.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_EDOM;
  }
  else if(n == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = 1.0 + a - x;
    result->err = 2.0 * GSL_DBL_EPSILON * (1.0 + fabs(a) + fabs(x));
    return GSL_SUCCESS;
  }
  else if(x < 0.0 || n < 5) {
    /* The explicit polynomial is always safe for x < 0
     * since all the terms are positive. Note that this
     * also catches overflows correctly.
     */
    return laguerre_n_cp(n, a, x, result);
  }
  else if(n > 1.0e+07 && x < 2.0*(a+1.0)+4.0*n) {
    return laguerre_large_n(n, a, x, result);
  }
  else {
    gsl_sf_result lg2;
    int stat_lg2 = gsl_sf_laguerre_2_impl(a, x, &lg2);
    double Lkm1 = 1.0 + a - x;
    double Lk   = lg2.val;
    double Lkp1;
    int k;

    for(k=2; k<n; k++) {
      Lkp1 = (-(k+a)*Lkm1 + (2.0*k+a+1.0-x)*Lk)/(k+1.0);
      Lkm1 = Lk;
      Lk   = Lkp1;
    }
    result->val = Lk;
    result->err = GSL_DBL_EPSILON * n * fabs(result->val);
    return stat_lg2;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_laguerre_1_e(double a, double x, gsl_sf_result * result)
{
  int status = gsl_sf_laguerre_1_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_laguerre_1_e", status);
  }
  return status;
}

int gsl_sf_laguerre_2_e(double a, double x, gsl_sf_result * result)
{
  int status = gsl_sf_laguerre_2_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_laguerre_2_e", status);
  }
  return status;
}

int gsl_sf_laguerre_3_e(double a, double x, gsl_sf_result * result)
{
  int status = gsl_sf_laguerre_3_impl(a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_laguerre_3_e", status);
  }
  return status;
}

int gsl_sf_laguerre_n_e(int n, double a, double x, gsl_sf_result * result)
{
  int status = gsl_sf_laguerre_n_impl(n, a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_laguerre_n_e", status);
  }
  return status;
}
