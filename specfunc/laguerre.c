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
laguerre_large_n(const int n, const double alpha, const double x, double * result)
{
  const double a = -n;
  const double b = alpha + 1.0;
  const double eta    = 2.0*b - 4.0*a;
  const double cos2th = x/eta;
  const double sin2th = 1.0 - cos2th;
  const double th = acos(sqrt(cos2th));
  const double pre_h  = 0.25*M_PI*M_PI*eta*eta*cos2th*sin2th;
  double ser;
  double lnpre;
  double lg_b;
  double lnfact;
  gsl_sf_lngamma_impl(b+n, &lg_b);
  gsl_sf_lnfact_impl(n, &lnfact);
  lnpre = lg_b - lnfact + 0.5*x + 0.5*(1.0-b)*log(0.25*x*eta) - 0.25*log(pre_h);
  ser = sin(a*M_PI) + sin(0.25*eta*(2.0*th - sin(2.0*th)) + 0.25*M_PI);
  return gsl_sf_exp_mult_impl(lnpre, ser, result);
}


/* Evaluate polynomial based on confluent hypergeometric representation.
 *
 * L^a_n(x) = (a+1)_n / n! 1F1(-n,a+1,x)
 *
 * assumes n > 0
 */
static
int
laguerre_n_cp(const int n, const double a, const double x, double * result)
{
  double lnfact;
  double lnpoch;
  int stat_f = gsl_sf_lnfact_impl(n, &lnfact);
  int stat_p = gsl_sf_lnpoch_impl(a+1.0, n, &lnpoch);
  double pre_correct;
  double lnpre;  double poly_1F1 = 1.0;
  int k;

/*
  if(fabs(lnpoch - lnfact) < 0.1*fabs(lnpoch + lnfact)) {
    lnpre = (lnpoch - M_LN2) - lnfact;
    pre_correct = 2.0;
  }
  else {
    lnpre = lnpoch - lnfact;
    pre_correct = 1.0;
  }
*/

  for(k=n-1; k>=0; k--) {
    double t = (-n+k)/(a+1.0+k) * (x/(k+1));
    double r = t + 1.0/poly_1F1;
    if(r > 0.9*GSL_DBL_MAX/poly_1F1) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    else {
      /* Collect the Horner terms. */
      poly_1F1 = 1.0 + t * poly_1F1;
    }
  }

  if(poly_1F1 == 0.0) {
    *result = 0.0;
    return GSL_SUCCESS;
  }
  else {
    int stat_e = gsl_sf_exp_mult_impl(lnpre, poly_1F1, result);
    *result *= pre_correct;
    return GSL_ERROR_SELECT_3(stat_e, stat_f, stat_p);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_laguerre_n_impl(const int n, const double a, const double x, double * result)
{
  if(n < 0 || a <= -1.0) {
    return GSL_EDOM;
  }
  else if(n == 0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    *result = 1.0 + a - x;
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
    double Lkm1 = 1.0 + a - x;
    double Lk   = gsl_sf_laguerre_2(a, x);
    double Lkp1;
    int k;

    for(k=2; k<n; k++) {
      Lkp1 = (-(k+a)*Lkm1 + (2.0*k+a+1.0-x)*Lk)/(k+1.0);
      Lkm1 = Lk;
      Lk   = Lkp1;
    }
    *result = Lk;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Error Handling Versions *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_laguerre_n_e(int n, double a, double x, double * result)
{
  int status = gsl_sf_laguerre_n_impl(n, a, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_laguerre_n_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*/

double gsl_sf_laguerre_1(const double a, const double x)
{
  return 1.0 + a - x;
}

double gsl_sf_laguerre_2(const double a, const double x)
{
  double c0 = 0.5 * (2.0+a)*(1.0+a);
  double c1 = -(2.0+a);
  double c2 = -0.5/(2.0+a);
  return c0 + c1*x*(1.0 + c2*x);
}

double gsl_sf_laguerre_3(const double a, const double x)
{
  double c0 = (3.0+a)*(2.0+a)*(1.0+a) / 6.0;
  double c1 = -c0 * 3.0 / (1.+a);
  double c2 = -1.0/(2.0+a);
  double c3 = -1.0/(3.0*(3.0+a));
  return c0 + c1*x*(1.0 + c2*x*(1.0 + c3*x));
}

double gsl_sf_laguerre_4(const double a, const double x)
{
  double c0 = (4.0+a)*(3.0+a)*(2.0+a)*(1.0+a) / 24.0;
  double c1 = -c0 * 4.0 / (1.0+a);
  double c2 = -3.0/(2.0*(2.0+a));
  double c3 = -2.0/(3.0*(3.0+a));
  double c4 = -1.0/(4.0*(4.0+a));
  return c0 + c1*x*(1.0 + c2*x*(1.0 + c3*x*(1.0 + c4*x)));
}

double gsl_sf_laguerre_5(const double a, const double x)
{
  double c0 = (5.0+a)*(4.0+a)*(3.0+a)*(2.0+a)*(1.0+a) / 120.0;
  double c1 = -c0 * 5.0 / (1.0+a);
  double c2 = -2.0/(2.0+a);
  double c3 = -1.0/(3.0+a);
  double c4 = -1.0/(2.0*(4.0+a));
  double c5 = -1.0/(5.0*(5.0+a));
  return c0 + c1*x*(1.0 + c2*x*(1.0 + c3*x*(1.0 + c4*x*(1.0 + c5*x))));
}

double gsl_sf_laguerre_n(int n, double a, double x)
{
  double y;
  int status = gsl_sf_laguerre_n_impl(n, a, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_laguerre_n", status);
  }
  return y;
}
