/* Author:  G. Jungman
 * RCS:     $Id$
 */
/* Miscellaneous implementations of use
 * for evaluation of hypergeometric functions.
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "hyperg.h"

#define locMAX(a,b)  ((a) > (b) ? (a) : (b))
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
