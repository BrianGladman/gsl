/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include "gsl_sf_hyperg.h"


static int hyperg_1F1_series(const double a, const double b, const double x,
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
  
  while(abs_del > GSL_MACH_EPS) {
    double u, abs_u;

    if(bn == 0.0) {
      return GSL_EDOM;
    }
    if(an == 0.0 || n > 200.0) {
      max_abs_del *= GSL_MACH_EPS;
      *prec   = fabs(GSL_MACH_EPS * n + max_abs_del);
      *result = sum;
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

    max_abs_del = locMAX(fabs(del), max_abs_del);

    an += 1.0;
    bn += 1.0;
    n  += 1.0;
  }
  
  max_abs_del *= GSL_MACH_EPS;
  *prec   = fabs(GSL_MACH_EPS * n + max_abs_del);
  *result = sum;
  return GSL_SUCCESS;
}

/* [Carlson, p.109] says the error in truncating this asymptotic series
 * is less than the absolute value of the first neglected term
 */
static int hyperg_2F0_series(const double a, const double b, const double x,
                             double * result, double * prec
			     )
{
  double an = a;
  double bn = b;
}

static int hyperg_1F1_asymp(const double a, const double b, const double x,
                            double * result, double * prec
			    )
{
  double pre, ln_pre, F;

  if(x == 0.0) {
    *prec   = 1.0;
    *result = 0.0; /* FIXME: ?? */
    return GSL_ELOSS;
  }
  else if(x < 0.0) {
    ln_pre = gsl_sf_lngamma(b) - a*log(-x) - gsl_sf_lngamma(b-a);
  }
  else {
    ln_pre = gsl_sf_lngamma(b) + x + (a-b)*log(x) - gsl_sf_lngamma(a);
  }
  
  
}
