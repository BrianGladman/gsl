/* integration/qagil.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "integration.h"

/* Evaluate an integral over an infinite range using the transformation

   integrate(f(x),-Inf,b) = integrate(f(b-(1-t)/t)/t^2,0,1)

   */

struct params { double b ; gsl_function * f ; } ;

static double transform (double t, void *params);

int
gsl_integration_qagil (gsl_function * f,
		       double b,
		       double epsabs, double epsrel, size_t limit,
		       gsl_integration_workspace * workspace,
		       double *result, double *abserr)
{
  int status;

  gsl_function f_transform;
  struct params transform_params  ;

  transform_params.b = b ;
  transform_params.f = f ;

  f_transform.function = &transform;
  f_transform.params = &transform_params;

  status = gsl_integration_qags_impl (&f_transform, 0.0, 1.0, 
				      epsabs, epsrel, limit,
				      workspace,
				      result, abserr,
				      &gsl_integration_qk15);

  return status;
}

static double 
transform (double t, void *params)
{
  struct params *p = (struct params *) params;
  double b = p->b;
  gsl_function * f = p->f;
  double x = b - (1 - t) / t;
  double y = GSL_FN_EVAL (f, x);
  return (y / t) / t;
}
