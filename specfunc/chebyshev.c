/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdlib.h>
#include <math.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_chebyshev.h"

int gsl_sf_cheb_calc_impl(struct gsl_sf_ChebSeries *, double (*)(double));


/*-*-*-*-*-*-*-*-*-*-*-* Allocators *-*-*-*-*-*-*-*-*-*-*-*/

struct gsl_sf_ChebSeries * gsl_sf_cheb_new(double (*func)(double),
    	    	    	    	    	   double a, double b,
			      	    	   int order)
{
  if(order < 0) {
    GSL_ERROR_RETURN("gsl_sf_cheb_new: order < 0", GSL_EDOM, 0);
  }
  else if(a >= b) {
    GSL_ERROR_RETURN("gsl_sf_cheb_new: null interval, a>=b", GSL_EDOM, 0);
  }
  else {
    struct gsl_sf_ChebSeries * cs = (struct gsl_sf_ChebSeries *)
      malloc(sizeof(struct gsl_sf_ChebSeries));
  
    if(cs == 0) {
      GSL_ERROR_RETURN("gsl_sf_cheb_new: out of memory", GSL_ENOMEM, 0);
    }

    cs->order = order;
    cs->a = a;
    cs->b = b;
    cs->c = malloc((order+1) * sizeof(double));
    if(cs->c == 0) {
      GSL_ERROR_RETURN("gsl_sf_cheb_new: out of memory", GSL_ENOMEM, 0);
    }
  
    if(gsl_sf_cheb_calc_impl(cs, func) != GSL_SUCCESS) {
      free(cs);
      return 0;
    }
    else {
      return cs;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_cheb_calc_impl(struct gsl_sf_ChebSeries * cs, double (*func)(double))
{
  if(cs == 0) {
    return GSL_EFAILED;
  }
  else {
    int k, j;
    double bma = 0.5 * (cs->b - cs->a);
    double bpa = 0.5 * (cs->b + cs->a);
    double fac = 2./(cs->order +1.);
    double * f = malloc((cs->order+1) * sizeof(double));

    if(f == 0) {
      return GSL_ENOMEM;
    }
  
    for(k = 0; k<=cs->order; k++) {
      double y = cos(M_PI * (k+0.5)/(cs->order+1));
      f[k] = func(y*bma + bpa);
    }

    for(j = 0; j<=cs->order; j++) {
      double sum = 0.0;
      for(k = 0; k<=cs->order; k++) sum += f[k]*cos(M_PI * j*(k+0.5)/(cs->order+1));
      cs->c[j] = fac * sum;
    }

    free(f);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_cheb_calc_e(struct gsl_sf_ChebSeries * cs, double (*func)(double))
{
  int status = gsl_sf_cheb_calc_impl(cs, func);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_cheb_calc_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/


double gsl_sf_cheb_eval_n(double x, int n, const struct gsl_sf_ChebSeries * cs)
{
  int j;
  double d  = 0.;
  double dd = 0.;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2. * y;

  int eval_order = (n < cs->order ? n : cs->order);
  
  for(j = eval_order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }
  return y*d - dd + 0.5 * cs->c[0];
}


double gsl_sf_cheb_eval(double x, const struct gsl_sf_ChebSeries * cs)
{
  int j;
  double d  = 0.;
  double dd = 0.;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2. * y;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }
  return y*d - dd + 0.5 * cs->c[0];
}


void gsl_sf_cheb_free(struct gsl_sf_ChebSeries * cs)
{
  if(cs != 0) {
    if(cs->c != 0) free(cs->c);
    free(cs);
  }
}
