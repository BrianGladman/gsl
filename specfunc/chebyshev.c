/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdlib.h>
#include <math.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include "gsl_sf_chebyshev.h"

#define gslMAX(a,b) ((a) > (b) ? (a) : (b))
#define gslMIN(a,b) ((a) < (b) ? (a) : (b))


int gsl_sf_cheb_calc_impl(struct gsl_sf_cheb_series *, double (*)(double));


/*-*-*-*-*-*-*-*-*-*-*-* Allocators *-*-*-*-*-*-*-*-*-*-*-*/

struct gsl_sf_cheb_series * gsl_sf_cheb_new(double (*func)(double),
    	    	    	    	    	   const double a, const double b,
			      	    	   const int order)
{
  int status;

  if(order < 0) {
    GSL_ERROR_RETURN("gsl_sf_cheb_new: order < 0", GSL_EDOM, 0);
  }
  else if(a >= b) {
    GSL_ERROR_RETURN("gsl_sf_cheb_new: null interval, a>=b", GSL_EDOM, 0);
  }
  else if(fabs(b-a) < 1.e+6 * GSL_MACH_EPS) {
    GSL_ERROR_RETURN("gsl_sf_cheb_new: interval close to null", GSL_EFAILED, 0);
  }
  else {
    struct gsl_sf_cheb_series * cs = (struct gsl_sf_cheb_series *)
      malloc(sizeof(struct gsl_sf_cheb_series));
  
    if(cs == 0) {
      GSL_ERROR_RETURN("gsl_sf_cheb_new: out of memory", GSL_ENOMEM, 0);
    }

    cs->cp = (double *)0;
    cs->ci = (double *)0;

    cs->order = order;
    cs->a = a;
    cs->b = b;
    cs->c = malloc((order+1) * sizeof(double));
    if(cs->c == 0) {
      GSL_ERROR_RETURN("gsl_sf_cheb_new: out of memory", GSL_ENOMEM, 0);
    }
  
    status = gsl_sf_cheb_calc_impl(cs, func);
    if(status != GSL_SUCCESS) {
      free(cs);
      GSL_ERROR_RETURN("gsl_sf_cheb_new: calc failed", status, 0);
    }
    else {
      return cs;
    }
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_cheb_calc_impl(struct gsl_sf_cheb_series * cs, double (*func)(double))
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

int gsl_sf_cheb_calc_deriv_impl(struct gsl_sf_cheb_series * cs)
{
  int j;
  int n = cs->order + 1;
  double con = 2. / (cs->b - cs->a);
  
  if(cs->cp != 0) free(cs->cp);
  cs->cp = (double *)malloc(n * sizeof(double));
  if(cs->cp == 0) return GSL_ENOMEM;

  cs->cp[n-1] = 0.0;
  
  if(n > 1) {
    cs->cp[n-2] = 2*(n-1) * cs->c[n-1];
    for(j = n-3; j>=0; j--) cs->cp[j] = cs->cp[j+2] + 2*(j+1) * cs->c[j+1];
    for(j = 0  ; j<n ; j++) cs->cp[j] *= con;
  }
  return GSL_SUCCESS;
}

int gsl_sf_cheb_calc_integ_impl(struct gsl_sf_cheb_series * cs)
{
  int n = cs->order + 1;
  double con = 0.25 * (cs->b - cs->a);
  
  if(cs->ci != 0) free(cs->ci);
  cs->ci = (double *)malloc(n * sizeof(double));
  if(cs->ci == 0) return GSL_ENOMEM;

  if(n == 1) {
    cs->ci[0] = 0.;
  }
  else if(n == 2) {
    cs->ci[1] = con * cs->c[0];
    cs->ci[0] = 2. * cs->ci[1];
  }
  else {
    double sum = 0.0;
    double fac = 1.0;
    int j;
    for(j=1; j<=n-2; j++) {
      cs->ci[j] = con * (cs->c[j-1] - cs->c[j+1])/j;
      sum += fac * cs->ci[j];
      fac = -fac;
    }
    cs->ci[n-1] = con * cs->c[n-2]/(n-1);
    sum += fac * cs->ci[n-1];
    cs->ci[0] = 2. * sum;
  }
  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_cheb_calc_e(struct gsl_sf_cheb_series * cs, double (*func)(double))
{
  int status = gsl_sf_cheb_calc_impl(cs, func);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_cheb_calc_e", status);
  }
  return status;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*-*/


double gsl_sf_cheb_eval_n(const double x, const int n, const struct gsl_sf_cheb_series * cs)
{
  int j;
  double d  = 0.;
  double dd = 0.;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2. * y;

  int eval_order = gslMIN(n, cs->order);
  eval_order     = gslMAX(eval_order, 0);
  
  for(j = eval_order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }
  return y*d - dd + 0.5 * cs->c[0];
}

double gsl_sf_cheb_eval(const double x, const struct gsl_sf_cheb_series * cs)
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

double gsl_sf_cheb_eval_deriv(const double x, struct gsl_sf_cheb_series * cs)
{
  int j;
  double d  = 0.;
  double dd = 0.;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2. * y;
  
  if(cs->cp == (double *)0) gsl_sf_cheb_calc_deriv_impl(cs);

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->cp[j];
    dd = temp;
  }
  return y*d - dd + 0.5 * cs->cp[0];
}

double gsl_sf_cheb_eval_integ(const double x, struct gsl_sf_cheb_series * cs)
{
  int j;
  double d  = 0.;
  double dd = 0.;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2. * y;
  
  if(cs->ci == (double *)0) gsl_sf_cheb_calc_integ_impl(cs);

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->ci[j];
    dd = temp;
  }
  return y*d - dd + 0.5 * cs->ci[0];
}

void gsl_sf_cheb_free(struct gsl_sf_cheb_series * cs)
{
  if(cs != 0) {
    if(cs->c != 0) free(cs->c);
    free(cs);
  }
}
