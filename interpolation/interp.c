/* interpolation/interp.c
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
#include <gsl/gsl_errno.h>
#include "gsl_interp.h"


int
gsl_interp_eval_e (const gsl_interp * obj,
		      const double xa[], const double ya[], double x,
		      gsl_interp_accel * a, double *y)
{
  return obj->eval (obj, xa, ya, x, a, y);
}

double
gsl_interp_eval (const gsl_interp * obj,
		 const double xa[], const double ya[], double x,
		 gsl_interp_accel * a)
{
  double y;
  int status = obj->eval (obj, xa, ya, x, a, &y);
  if (status != GSL_SUCCESS)
    {
      GSL_WARNING ("gsl_interp_eval", status);
    }
  return y;
}


int
gsl_interp_eval_deriv_e (const gsl_interp * obj,
			    const double xa[], const double ya[], double x,
			    gsl_interp_accel * a,
			    double *dydx)
{
  return obj->eval_d (obj, xa, ya, x, a, dydx);
}

double
gsl_interp_eval_deriv (const gsl_interp * obj,
		       const double xa[], const double ya[], double x,
		       gsl_interp_accel * a)
{
  double dydx;
  int status = obj->eval_d (obj, xa, ya, x, a, &dydx);
  if (status != GSL_SUCCESS)
    {
      GSL_WARNING ("gsl_interp_eval_deriv", status);
    }
  return dydx;
}


int
gsl_interp_eval_deriv2_e (const gsl_interp * obj,
			     const double xa[], const double ya[], double x,
			     gsl_interp_accel * a,
			     double * d2)
{
  return obj->eval_d2 (obj, xa, ya, x, a, d2);
}

double
gsl_interp_eval_deriv2 (const gsl_interp * obj,
		        const double xa[], const double ya[], double x,
		        gsl_interp_accel * a)
{
  double d2;
  int status = obj->eval_d2 (obj, xa, ya, x, a, &d2);
  if (status != GSL_SUCCESS)
    {
      GSL_WARNING ("gsl_interp_eval_deriv2", status);
    }
  return d2;
}


int
gsl_interp_eval_integ_e (const gsl_interp * obj,
			    const double xa[], const double ya[],
                            double a, double b,
			    gsl_interp_accel * acc,
			    double * result)
{
  return obj->eval_i (obj, xa, ya, acc, a, b, result);
}


double
gsl_interp_eval_integ (const gsl_interp * obj,
		       const double xa[], const double ya[],
                       double a, double b,
		       gsl_interp_accel * acc)
{
  double result;
  int status = obj->eval_i (obj, xa, ya, acc, a, b, &result);
  if (status != GSL_SUCCESS)
    {
      GSL_WARNING ("gsl_interp_eval_integ", status);
    }
  return result;
}


void
gsl_interp_free (gsl_interp * obj)
{
  if (obj != 0)
    obj->free (obj);
}
