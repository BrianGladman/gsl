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
 * is less than the absolute value of the first neglected term.
 *
 * A termination argument is provided, so that the series will
 * be summed at most up to n=n_trunc. If n_trunc is set negative,
 * then the series is summed until it appears to start diverging.
 */
static
int
hyperg_2F0_series(const double a, const double b, const double x,
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


int
gsl_sf_hyperg_2F0_series_impl(const double a, const double b, const double x,
                              int n_trunc,
                              double * result
                              )
{
  double prec;
  return hyperg_2F0_series(a, b, x, n_trunc, result, &prec);
}


int
gsl_sf_hyperg_2F0_impl(const double a, const double b, const double x, double * result)
{
  if(x < 0.0) {
    /* Use "definition" 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)
     */
    double U;
    double pre = pow(-1.0/x, a);
    int stat_U = gsl_sf_hyperg_U_impl(a, 1.0+a-b, -1.0/x, &U);
    *result = pre * U;
    return stat_U;
  }
  else if(x == 0.0) {
    *result = 1.0;
    return GSL_SUCCESS;
  }
  else {
    /* Use asymptotic series. ??
     */
    /* return hyperg_2F0_series(a, b, x, -1, result, &prec); */
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
