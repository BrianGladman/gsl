/* specfunc/chebyshev.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "gsl_sf_chebyshev.h"


/*-*-*-*-*-*-*-*-*-*-*-* Allocators *-*-*-*-*-*-*-*-*-*-*-*/

gsl_sf_cheb_series * gsl_sf_cheb_new(double (*func)(double),
    	    	    	    	     const double a, const double b,
			      	     const int order)
{
  int status;

  if(order < 0) {
    GSL_ERROR_VAL("gsl_sf_cheb_new: order < 0", GSL_EDOM, 0);
  }
  else if(a >= b) {
    GSL_ERROR_VAL("gsl_sf_cheb_new: null interval, a>=b", GSL_EDOM, 0);
  }
#if 0
  else if(fabs(b-a) < 1.0e+6 * GSL_DBL_EPSILON) {
    /* FIXME: arbitrary nonsense */
    GSL_ERROR_VAL("gsl_sf_cheb_new: interval close to null", GSL_EFAILED, 0);
  }
#endif
  else {
    gsl_sf_cheb_series * cs = (gsl_sf_cheb_series *)
      malloc(sizeof(gsl_sf_cheb_series));
  
    if(cs == 0) {
      GSL_ERROR_VAL("gsl_sf_cheb_new: out of memory", GSL_ENOMEM, 0);
    }

    cs->cp = (double *)0;
    cs->ci = (double *)0;

    cs->order    = order;
    cs->order_sp = order;
    cs->a = a;
    cs->b = b;
    cs->c = (double *) malloc((order+1) * sizeof(double));
    if(cs->c == 0) {
      GSL_ERROR_VAL("gsl_sf_cheb_new: out of memory", GSL_ENOMEM, 0);
    }

    status = gsl_sf_cheb_calc_e(cs, func);
    if(status != GSL_SUCCESS) {
      free(cs);
      GSL_ERROR_VAL("gsl_sf_cheb_new: calc failed", status, 0);
    }
    else {
      return cs;
    }
  }
}


void gsl_sf_cheb_free(gsl_sf_cheb_series * cs)
{
  if(cs != 0) {
    if(cs->c != 0) free(cs->c);
    if(cs->cp != 0) free(cs->cp);
    if(cs->ci != 0) free(cs->ci);
    free(cs);
  }
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_cheb_calc_e(gsl_sf_cheb_series * cs, double (*func)(double))
{
  /* CHECK_POINTER(cs) */ 

  {
    int k, j;
    double bma = 0.5 * (cs->b - cs->a);
    double bpa = 0.5 * (cs->b + cs->a);
    double fac = 2.0/(cs->order +1.0);
    double * f = (double *) malloc((cs->order+1) * sizeof(double));

    if(f == 0) {
      GSL_ERROR ("error", GSL_ENOMEM);
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

int gsl_sf_cheb_calc_deriv_e(gsl_sf_cheb_series * cs)
{
  const int n = cs->order + 1;
  const double con = 2.0 / (cs->b - cs->a);
  int j;
  
  if(cs->cp != 0) free(cs->cp);
  cs->cp = (double *) malloc(n * sizeof(double));
  if(cs->cp == 0) GSL_ERROR ("error", GSL_ENOMEM);

  cs->cp[n-1] = 0.0;
  
  if(n > 1) {
    cs->cp[n-2] = 2*(n-1) * cs->c[n-1];
    for(j = n-3; j>=0; j--) cs->cp[j] = cs->cp[j+2] + 2*(j+1) * cs->c[j+1];
    for(j = 0  ; j<n ; j++) cs->cp[j] *= con;
  }
  return GSL_SUCCESS;
}

int gsl_sf_cheb_calc_integ_e(gsl_sf_cheb_series * cs)
{
  const int n = cs->order + 1;
  const double con = 0.25 * (cs->b - cs->a);
  
  if(cs->ci != 0) free(cs->ci);
  cs->ci = (double *) malloc(n * sizeof(double));
  if(cs->ci == 0) GSL_ERROR ("error", GSL_ENOMEM);

  if(n == 1) {
    cs->ci[0] = 0.;
  }
  else if(n == 2) {
    cs->ci[1] = con * cs->c[0];
    cs->ci[0] = 2.0 * cs->ci[1];
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
    cs->ci[0] = 2.0 * sum;
  }
  return GSL_SUCCESS;
}


int
gsl_sf_cheb_eval_n_e(const gsl_sf_cheb_series * cs,
                        const int n, const double x,
                        gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  int eval_order = GSL_MAX_INT(0, GSL_MIN_INT(n, cs->order));

  for(j = eval_order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  /* CHECK_POINTER(result) */

  {
    result->val = y*d - dd + 0.5 * cs->c[0];
    result->err = fabs(cs->c[eval_order]);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_cheb_eval_e(const gsl_sf_cheb_series * cs,
                      const double x,
                      gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  /* CHECK_POINTER(result) */

  {
    result->val = y*d - dd + 0.5 * cs->c[0];
    result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->c[cs->order]);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_cheb_eval_mode_e(const gsl_sf_cheb_series * cs,
                           const double x,
                           gsl_mode_t mode,
                           gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  int eval_order;

  if(GSL_MODE_PREC(mode) == GSL_PREC_DOUBLE)
    eval_order = cs->order;
  else
    eval_order = cs->order_sp;

  for(j = eval_order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  /* CHECK_POINTER(result) */

  {
    result->val = y*d - dd + 0.5 * cs->c[0];
    result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->c[eval_order]);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_cheb_eval_deriv_e(gsl_sf_cheb_series * cs, const double x,
                            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  if(cs->cp == (double *)0) gsl_sf_cheb_calc_deriv_e(cs);

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->cp[j];
    dd = temp;
  }

  /* CHECK_POINTER(result) */

  {
    result->val = y*d - dd + 0.5 * cs->cp[0];
    result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->cp[cs->order]);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_cheb_eval_integ_e(gsl_sf_cheb_series * cs, const double x,
                            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;
  
  if(cs->ci == (double *)0) gsl_sf_cheb_calc_integ_e(cs);

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->ci[j];
    dd = temp;
  }

  /* CHECK_POINTER(result) */

  {
    result->val = y*d - dd + 0.5 * cs->ci[0];
    result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->ci[cs->order]);
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_cheb_eval(const gsl_sf_cheb_series * cs, double x)
{
  EVAL_RESULT(gsl_sf_cheb_eval_e(cs, x, &result));
}

double gsl_sf_cheb_eval_n(const gsl_sf_cheb_series * cs, int order, double x)
{
  EVAL_RESULT(gsl_sf_cheb_eval_n_e(cs, order, x, &result));
}

double gsl_sf_cheb_eval_mode(const gsl_sf_cheb_series * cs,
                             const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_cheb_eval_n_e(cs, x, mode, &result));
}

double gsl_sf_cheb_eval_deriv(gsl_sf_cheb_series * cs, double x)
{
  EVAL_RESULT(gsl_sf_cheb_eval_deriv_e(cs, x, &result));
}

double gsl_sf_cheb_eval_integ(gsl_sf_cheb_series * cs, double x)
{
  EVAL_RESULT(gsl_sf_cheb_eval_integ_e(cs, x, &result));
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*/
