/* ode-initval/gear1.c
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
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


typedef  struct gsl_odeiv_step_gear1_struct  gsl_odeiv_step_gear1;

struct gsl_odeiv_step_gear1_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  double * k;
  double * y0;
};

static int  gear1_step(void *, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static void gear1_free(void *);



gsl_odeiv_step *
gsl_odeiv_step_gear1_new(void)
{
  gsl_odeiv_step_gear1 * s = (gsl_odeiv_step_gear1 *) malloc(sizeof(gsl_odeiv_step_gear1));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "gear1",
      gear1_step,
      0,
      gear1_free,
      0,
      0,
      2);

    s->k = 0;
    s->y0 = 0;
  }
  return (gsl_odeiv_step *) s;
}


static int
gear1_step(
  void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  const int iter_steps = 3;
  int status = 0;
  int nu;
  int i;
  size_t dim;

  gsl_odeiv_step_gear1 * my = (gsl_odeiv_step_gear1 *) self;

  if(sys->dimension <= 0) {
    return GSL_EINVAL;
  }

  /* (Re)Initialization may be necessary. */
  if(sys->dimension != my->parent.dimension) {
    if(my->k != 0) free(my->k);
    if(my->y0 != 0) free(my->y0);
    my->parent.dimension = sys->dimension;
    my->k  = (double *) malloc(sys->dimension * sizeof(double));
    my->y0 = (double *) malloc(sys->dimension * sizeof(double));
    if(my->k == 0 || my->y0 == 0) {
      if(my->k != 0) free(my->k);
      if(my->y0 != 0) free(my->y0);
      my->parent.dimension = 0;
      return GSL_ENOMEM;
    }
  }

  dim = my->parent.dimension;

  memcpy(my->y0, y, dim * sizeof(double));

  /* iterative solution */
  for(nu=0; nu<iter_steps; nu++) {
    status += ( GSL_ODEIV_FN_EVAL(sys, t + h, y, my->k) != 0 );
    for(i=0; i<dim; i++) {
      y[i] = my->y0[i] + h * my->k[i];
    }
  }

  /* fudge the error estimate */
  for(i=0; i<dim; i++) {
    yerr[i] = h * h * my->k[i];
  }

  if(dydt_out != 0) {
    memcpy(dydt_out, my->k, dim * sizeof(double));
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static void
gear1_free(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_gear1 * my = (gsl_odeiv_step_gear1 *) self;
    if(my->k != 0) free(my->k);
    if(my->y0 != 0) free(my->y0);
    free(self);
  }
}
