/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Miscellaneous implementations of use
 * for evaluation of hypergeometric functions.
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "hyperg.h"

#define locMAX(a,b)  ((a) > (b) ? (a) : (b))
#define locMIN(a,b)  ((a) < (b) ? (a) : (b))
#define locEPS       (1000.0*GSL_MACH_EPS)
#define SUM_LARGE    (1.0e-5*DBL_MAX)


int
gsl_sf_hyperg_1F1_series_impl(const double a, const double b, const double x,
                              double * result, double * prec
                              )
{
  double an  = a;
  double bn  = b;
  double n   = 1.0;
  double sum = 1.0;
  double del = 1.0;
  double abs_del = 1.0;
  double max_abs_del = 1.0;
  double err;
  
  while(abs_del/fabs(sum) > GSL_MACH_EPS) {
    double u, abs_u;

    if(bn == 0.0) {
      *result = 0.0;
      return GSL_EDOM;
    }
    if(an == 0.0 || n > 200.0) {
      max_abs_del *= GSL_MACH_EPS;
      err     = fabs(GSL_MACH_EPS * n + max_abs_del);
      *prec   = err/(err + fabs(sum));
      *result = sum;
      if(*prec > locEPS)
        return GSL_ELOSS;
      else
        return GSL_SUCCESS;
    }

    u = x * (an/(bn*n));
    abs_u = fabs(u);
    if(abs_u > 1.0 && max_abs_del > DBL_MAX/abs_u) {
      *prec   = 1.0;
      *result = sum;
      return GSL_ELOSS;
    }
    del *= u;
    sum += del;

    if(fabs(sum) > SUM_LARGE) {
      *prec = 1.0;
      *result = sum;
      return GSL_EOVRFLW;
    }

    abs_del = fabs(del);
    max_abs_del = locMAX(abs_del, max_abs_del);

    an += 1.0;
    bn += 1.0;
    n  += 1.0;
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


int
gsl_sf_hyperg_1F1_large_b_impl(const double a, const double b, const double x, double * result)
{
  if(fabs(x/b) < 1.0) {
    const double u = x/b;
    const double v = 1.0/(1.0-u);
    const double pre = pow(v,a);
    const double uv  = u*v;
    const double uv2 = uv*uv;
    const double t1  = a*(a+1.0)/(2.0*b)*uv2;
    const double t2a = a*(a+1.0)/(24.0*b*b)*uv2;
    const double t2b = 12.0 + 16.0*(a+2.0)*uv + 3.0*(a+2.0)*(a+3.0)*uv2;
    const double t2  = t2a*t2b;
    *result = pre * (1.0 - t1 + t2);
    return GSL_SUCCESS;
  }
  else {
    return GSL_EDOM;
  }
}


int
gsl_sf_hyperg_U_large_b_impl(const double a, const double b, const double x, double * result)
{
  double N   = floor(b);  /* b = N + eps */
  double eps = b - N;
  
  if(fabs(eps) < GSL_SQRT_MACH_EPS) {
    double lnpre;
    double M;
    if(b > 0.0) {
      double lg_bm1;
      double lg_a;
      gsl_sf_lngamma_impl(b-1.0, &lg_bm1);
      gsl_sf_lngamma_impl(a, &lg_a);
      lnpre = (1.0-b)*log(x) + x + lg_bm1 - lg_a;
      gsl_sf_hyperg_1F1_large_b_impl(1.0-a, 2.0-b, -x, &M);
    }
    else {
      double lg_1mb;
      double lg_1pamb;
      gsl_sf_lngamma_impl(1.0-b, &lg_1mb);
      gsl_sf_lngamma_impl(1.0+a-b, &lg_1pamb);
      lnpre = lg_1mb - lg_1pamb;
      gsl_sf_hyperg_1F1_large_b_impl(a, b, x, &M);
    }
    if(lnpre > GSL_LOG_DBL_MAX-10.0) {
      *result = 0.0;
      return GSL_EOVRFLW;
    }
    else {
      *result = exp(lnpre) * M;
      return GSL_SUCCESS;
    }
  }
  else {
    double lg_1mb,   sgn_1mb;
    double lg_1pamb, sgn_1pamb;
    double lg_bm1,   sgn_bm1;
    double lg_a,     sgn_a;
    double M1, M2;
    double lnpre1, lnpre2;
    gsl_sf_hyperg_1F1_large_b_impl(    a,     b, x, &M1);
    gsl_sf_hyperg_1F1_large_b_impl(1.0-a, 2.0-b, x, &M2);

    gsl_sf_lngamma_sgn_impl(1.0-b,   &lg_1mb,   &sgn_1mb);
    gsl_sf_lngamma_sgn_impl(1.0+a-b, &lg_1pamb, &sgn_1pamb);

    gsl_sf_lngamma_sgn_impl(b-1.0, &lg_bm1, &sgn_bm1);
    gsl_sf_lngamma_sgn_impl(a,     &lg_a,   &sgn_a);

    lnpre1 = lg_1mb - lg_1pamb;
    lnpre2 = lg_bm1 - lg_a - (1.0-b)*log(x) - x;

    if(lnpre1 > GSL_LOG_DBL_MAX-10.0 || lnpre2 > GSL_LOG_DBL_MAX-10.0) {
      *result = 0.0;
      return GSL_EOVRFLW;
    }
    else {
      *result = exp(lnpre1) * M1 + exp(lnpre2)*M2;
      return GSL_SUCCESS;
    }
  }
}



/* [Carlson, p.109] says the error in truncating this asymptotic series
 * is less than the absolute value of the first neglected term.
 *
 * A termination argument is provided, so that the series will
 * be summed at most up to n=n_trunc. If n_trunc is set negative,
 * then the series is summed until it appears to start diverging.
 */
int
gsl_sf_hyperg_2F0_series_impl(const double a, const double b, const double x,
                              int n_trunc,
                              double * result, double * prec
                              )
{
  double an = a;
  double bn = b;  
  double n   = 1.0;
  double sum = 1.0;
  double del = 1.0;
  double abs_del = 1.0;
  double max_abs_del = 1.0;
  double last_abs_del = 1.0;
  double err;
  
  while(abs_del/fabs(sum) > GSL_MACH_EPS && n < 200.0) {

    double u = an * (bn/n * x);
    double abs_u = fabs(u);

    if(abs_u > 1.0 && (max_abs_del > DBL_MAX/abs_u)) {
      *prec = 1.0;
      *result = sum;
      return GSL_EOVRFLW;
    }

    del *= u;
    sum += del;

    abs_del = fabs(del);

    if(abs_del > last_abs_del) break; /* series is probably starting to grow */

    last_abs_del = abs_del;
    max_abs_del  = locMAX(abs_del, max_abs_del);

    an += 1.0;
    bn += 1.0;
    n  += 1.0;
    
    if(an == 0.0 || bn == 0.0) break;        /* series terminated */
    
    if(n_trunc >= 0 && n >= n_trunc) break;  /* reached requested timeout */
  }

  err     = GSL_MACH_EPS * n + abs_del;
  *prec   = err/(err + fabs(sum));
  *result = sum;
  if(*prec > locEPS)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}
