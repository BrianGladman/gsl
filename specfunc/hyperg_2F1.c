/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_gamma.h"
#include "gsl_sf_hyperg.h"

#define locMAX(a, b)    ((a) > (b) ? (a) : (b))
#define locEPS          (1000.0*GSL_MACH_EPS)


static int hyperg_2F1_series(const double a, const double b, const double c,
                             const double x, 
			     double * result, double * prec
			     )
{
  double sum = 1.0;
  double del = 1.0;
  double delmax = 0.0;
  double k = 0.0;
  int i = 0;

  if(fabs(c) < locEPS) {
    *prec   = 1.0;
    *result = 0.0; /* FIXME: ?? */
    return GSL_ELOSS;
  }

  do {
    if(++i > 20000) {
      *prec   = 1.0;
      *result = sum;
      return GSL_ELOSS;
    }
    del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  /* Gauss series */
    sum += del;
    delmax = locMAX(fabs(del), delmax);
    k += 1.0;
  } while(fabs(del/sum) > GSL_MACH_EPS);

  *prec   = (GSL_MACH_EPS*delmax)/fabs(sum) + i*GSL_MACH_EPS;
  *result = sum;
  return GSL_SUCCESS;
}

static int pow_omx(const double x, const double p, double * result)
{
  double ln_omx;
  double ln_result;
  if(fabs(x) < GSL_ROOT3_MACH_EPS) {
    ln_omx = -x*(1.0 + 0.5*x + x*x/3.0);
  }
  else {
    ln_omx = log(1.0-x);
  }
  ln_result = p * ln_omx;
  if(ln_result > GSL_LOG_DBL_MAX) {
    *result = 0.0; /* FIXME: should be Inf */
    return GSL_EOVRFLW;
  }
  else if(ln_result < GSL_LOG_DBL_MIN) {
    *result = 0.0;
    return GSL_EUNDRFLW;
  }
  else {
    *result = exp(ln_result);
    return GSL_SUCCESS;
  }
}

int gsl_sf_hyperg_2F1_impl(const double a, const double b, const double c,
                           const double x,
			   double * result
			   )
{
  const double ax = fabs(x);
  const double d  = c - a - b;
  const int a_neg_integer = ( a < 0.0  &&  fabs(a - rint(a)) < locEPS );
  const int b_neg_integer = ( b < 0.0  &&  fabs(b - rint(b)) < locEPS );
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rint(c)) < locEPS );
  const int d_integer     = ( fabs(d - rint(d)) < locEPS );
  
  if(ax >= 1.0) return GSL_EDOM;

  if(c_neg_integer) {
    if(! (a_neg_integer && a > c)) return GSL_EDOM;
    if(! (b_neg_integer && b > c)) return GSL_EDOM;
  }
  
  if(fabs(c-b) < locEPS || fabs(c-a) < locEPS) {
    return pow_omx(x, c-a-b, result);  /* (1-x)^(c-a-b) */
  }
  
  if(x < -0.5 && b > 0.0) {
    double F, F_prec, p;
    int status_F = hyperg_2F1_series(a, c-b, c, -x/(1.0-x), &F, &F_prec);
    int status_p = pow_omx(x, -a, &p);
    if(status_p == GSL_SUCCESS && status_F == GSL_SUCCESS) {
      *result = p * F;
      return GSL_SUCCESS; 
    }
    else if(status_p == GSL_EUNDRFLW) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else if(status_p == GSL_EOVRFLW) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    else if(F_prec > locEPS) {
      *result = p * F;
      return GSL_ELOSS;
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
  
  if(x < -0.5 && b <= 0.0) {
    double F, F_prec, p;
    int status_F = hyperg_2F1_series(c-a, b, c, -x/(1.0-x), &F, &F_prec);
    int status_p = pow_omx(x, -b, &p);
    if(status_p == GSL_SUCCESS && status_F == GSL_SUCCESS) {
      *result = p * F;
      return GSL_SUCCESS; 
    }
    else if(status_p == GSL_EUNDRFLW) {
      *result = 0.0;
      return GSL_EUNDRFLW;
    }
    else if(status_p == GSL_EOVRFLW) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    else if(F_prec > locEPS) {
      *result = p * F;
      return GSL_ELOSS;
    }
    else {
      *result = 0.0;
      return GSL_EFAILED;
    }
  }
  
  if(x > 0.75 && !d_integer) {
    double pre1, pre2, F1, F2;
    double status_F1, status_F2;
    double prec_F1, prec_F2;
    double lng_c   = gsl_sf_lngamma(c);
    double ln_pre1 = lng_c + gsl_sf_ln_gamma( d) - gsl_sf_lngamma(c-a) - gsl_sf_lngamma(c-b);
    double ln_pre2 = lng_c + gsl_sf_ln_gamma(-d) - gsl_sf_lngamma(a)   - gsl_sf_lngamma(b) + d*log(1.0-x);
    if(ln_pre1 > GSL_LOG_DBL_MAX || ln_pre2 > GSL_LOG_DBL_MAX) {
      *result = 0.0; /* FIXME: should be Inf */
      return GSL_EOVRFLW;
    }
    pre1 = exp(ln_pre1);
    pre2 = exp(ln_pre2);
    status_F1 = hyperg_2F1_series(  a,   b, 1.0-d, 1.0-x, &F1, &prec_F1);
    status_F2 = hyperg_2F1_series(c-a, c-b, 1.0+d, 1.0-x, &F2, &prec_F2);
    *result = pre1*F1 + pre2*F2;
    if(prec_F1 > locEPS || prec_F2 > locEPS) {
      return GSL_ELOSS;
    }
    else {
      return GSL_SUCCESS;
    }
  }
  
  if(x > 0.75 && d_integer) {
    double ad = fabs(d);
    double ln_pre1, ln_pre2, pre1, pre2;
    double F1, F2;
    double d1, d2;
    if(d >= 0.0) {
      d1 = d;
      d2 = 0.0;
    }
    else {
      d1 = 0.0;
      d2 = d;
    }
    ln_pre1 = gsl_sf_lngamma(ad) + gsl_sf_lngamma(c) + d2*log(1.0-x) - gsl_sf_lngamma(a+d1) - gsl_sf_lngamma(b+d1);
    ln_pre2 = gsl_sf_lngamma(c) + d1*log(1.0-x) - gsl_sf_lngamma(a+d2) - gsl_sf_lngamma(b+d2);
    if(ln_pre1 > GSL_LOG_DBL_MAX || ln_pre2 > GSL_LOG_DBL_MAX) {
      *result = 0.0; /* FIXME: should be Inf (?) */
      return GSL_EOVRFLW;
    }
    pre1 = exp(ln_pre1);
    pre2 = exp(ln_pre2);
    
    /* FIXME: F1= ...  F2= ... */
  }
  
  /* default case */
  {
    double prec;
    int status = hyperg_2F1_series(a, b, c, x, result, &prec);
    if(prec > locEPS) {
      return GSL_ELOSS;
    }
    else {
      return GSL_SUCCESS;
    }
  }
}
