/* ode-initval/odeiv.c
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
#include <gsl/gsl_errno.h>
#include "gsl_odeiv.h"


const char *
gsl_odeiv_step_name(const gsl_odeiv_step * s)
{
  return s->_name;
}


int
gsl_odeiv_step_apply(
  gsl_odeiv_step * step,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * dydt)
{
  return step->_step(step, t, h, y, yerr, dydt_in, dydt_out, dydt);
}


int
gsl_odeiv_step_reset(gsl_odeiv_step * s)
{
  if(s->_reset != 0) {
    return s->_reset(s);
  }
  else {
    return GSL_SUCCESS;
  }
}


void
gsl_odeiv_step_free(gsl_odeiv_step * s)
{
  if(s != 0) {
    if(s->_free != 0)  s->_free(s);
  }
}
