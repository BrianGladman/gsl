/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_hyperg.h"

#define locEPS         (1000.0 * GSL_MACH_EPS)
#define locMAX(a, b)   ((a) > (b) ? (a) : (b))


/* [Carlson, p.109] says the error in truncating this asymptotic series
 * is less than the absolute value of the first neglected term
 */
static
int
hyperg_2F0_series(const double a, const double b, const double x,
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

    abs_del      = fabs(del);

    if(abs_del > last_abs_del) break; /* series is probably starting to grow */

    last_abs_del = abs_del;
    max_abs_del  = locMAX(abs_del, max_abs_del);

    an += 1.0;
    bn += 1.0;
    n  += 1.0;
    
    if(an == 0.0 || bn == 0.0) break;  /* series terminated */
  }

  err     = GSL_MACH_EPS * n + abs_del;
  *prec   = err/(err + fabs(sum));
  *result = sum;
  if(*prec > locEPS)
    return GSL_ELOSS;
  else
    return GSL_SUCCESS;
}


int
gsl_sf_hyperg_2F0_impl(const double a, const double b, const double x,
                       double * result
                       )
{
  if(fabs(x) < 0.1 && fabs(a) < 1.0 && fabs(b) < 1.0) {
    double prec;
    return hyperg_2F0_series(a, b, x, result, &prec);
  }
  else {
    /* FIXME: can we do something here? */
    *result = 0.0;
    return GSL_EDOM;
  }
}


int
gsl_sf_hyperg_2F0_e(const double a, const double b, const double x, double * result)
{
  int status = gsl_sf_hyperg_2F0_impl(a, b, x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_hyperg_2F0_e", status);
  }
  return status;
}

double
gsl_sf_hyperg_2F0(const double a, const double b, const double x)
{
  double y;
  int status = gsl_sf_hyperg_2F0_impl(a, b, x, &y);
  if(status != GSL_SUCCESS) {
    GSL_WARNING("gsl_sf_hyperg_2F0", status);
  }
  return y;
}


int
test_hyperg2F0_stuff(void)
{
  double x = -0.1;
  double a = 1.0;
  double b = 1.0;
  double p;
  double r;
  
  hyperg_2F0_series(a, b, x, &r, &p);
  
  printf("%9.6g  %9.6g  %9.6g    %20.16g  %20.16g\n", a, b, x, r, p);
  
  return 0;
}
