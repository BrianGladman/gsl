/* ode-initval/rk4.c
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
#include <gsl/gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


typedef  struct gsl_odeiv_step_rk4_struct  gsl_odeiv_step_rk4;

struct gsl_odeiv_step_rk4_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  double * work;  /* generic work space */
};


static int  rk4_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static void rk4_free(void * self);


gsl_odeiv_step *
gsl_odeiv_step_rk4_new(void)
{
  gsl_odeiv_step_rk4 * s = (gsl_odeiv_step_rk4 *) malloc(sizeof(gsl_odeiv_step_rk4));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "rk4",
      rk4_step,
      0,
      rk4_free,
      0,
      0,
      4);
    s->work = 0;
  }
  return (gsl_odeiv_step *) s;
}


static int
rk4_step(void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  int i;
  int status = 0;
  size_t dim;

  double * k;
  double * y0;
  double * ytmp;
 
  gsl_odeiv_step_rk4 * my = (gsl_odeiv_step_rk4 *) self;

  if(sys->dimension <= 0) {
    return GSL_EINVAL;
  }  

  if(sys->dimension != my->parent.dimension) {
    if(my->work != 0) free(my->work);
    my->parent.dimension = sys->dimension;
    my->work = (double *) malloc(3 * sys->dimension * sizeof(double));
    if(my->work == 0) {
      my->parent.dimension = 0;
      return GSL_ENOMEM;
    }
  }

  dim = my->parent.dimension;

  /* divide up the workspace */  
  k    = my->work;
  y0   = my->work + dim;
  ytmp = my->work + 2*dim;

  /* Copy the starting value. We will write over
   * the y[] vector, using it for scratch and
   * then filling it with the final result.
   */
  memcpy(y0, y, dim * sizeof(double));

  /* k1 step */
  if(dydt_in != 0) {
    memcpy(k, dydt_in, dim * sizeof(double));
  }
  else {
    status += ( GSL_ODEIV_FN_EVAL(sys, t, y0, k) != 0 );
  }
  for(i=0; i<dim; i++) {
    y[i] = h/6.0 * k[i]; /* use y[] to store delta_y */
    ytmp[i] = y0[i] + 0.5*h * k[i];
  }

  /* k2 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + 0.5*h, ytmp, k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/3.0 * k[i];
    ytmp[i] = y0[i] + 0.5*h * k[i];
  }  

  /* k3 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + 0.5*h, ytmp, k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/3.0 * k[i];
    ytmp[i] = y0[i] + h * k[i];
  } 

  /* k4 step, error estimate, and final sum */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + h, ytmp, k) != 0 );
  for(i=0; i<dim; i++) {
    y[i]   += h/6.0 * k[i];
    yerr[i] = h * y[i];
    y[i]   += y0[i];
    if(dydt_out != 0) dydt_out[i] = k[i];
  } 

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static void
rk4_free(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_rk4 * my = (gsl_odeiv_step_rk4 *) self;
    if(my->work != 0) free(my->work);
    free(self);
  }
}
