/* ode-initval/odeiv_util.c
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
#include <string.h>
#include "odeiv_util.h"

/*
gsl_odeiv_step *
gsl_odeiv_step_new(
  const char * name,
  unsigned int dim,
  unsigned int ord,
  size_t state_size,
  size_t work_size)
{
  gsl_odeiv_step * s;

  if(dim == 0 || name == 0) return 0;

  s = (gsl_odeiv_step *) malloc(sizeof(gsl_odeiv_step));

  if(s != 0) {
    s->_name = (char *) malloc(strlen(name));

    s->dimension = dim;
    s->order = ord;

    s->_step  = 0;
    s->_reset = 0;
    s->_free  = 0;
    s->_state = 0;
    s->_work  = 0;

    if(state_size > 0) s->_state = malloc(state_size);
    if(work_size > 0) s->_work = malloc(work_size);

    if((s->_state == 0 && state_size > 0) || (s->_work == 0 && work_size > 0) || s->_name == 0) {
      if(s->_name != 0) free(s->_name);
      if(s->_state != 0) free(s->_state);
      if(s->_work != 0) free(s->_work);
      free(s);
      return 0;
    }

    strcpy(s->_name, name);
  }

  return s;
}
*/

void gsl_odeiv_step_construct(
  gsl_odeiv_step * s,
  const char * _name,
  int  (*_step)  (void *, double, double, double *, double *, const double *, double *, const gsl_odeiv_system *),
  int  (*_reset) (void *),
  void (*_free)  (void *),
  int can_use_dydt,
  size_t dimension,
  unsigned int order
  )
{
  s->_name = _name;
  s->_step = _step;
  s->_reset = _reset;
  s->_free = _free;
  s->can_use_dydt = can_use_dydt;
  s->dimension = dimension;
  s->order = order;
}
