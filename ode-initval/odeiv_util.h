/* ode-initval/odeiv_util.h
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
#ifndef ODEIV_UTIL_H
#define ODEIV_UTIL_H

#include <stdlib.h>
#include "gsl_odeiv.h"


/* Low-level allocator for gsl_odeiv_step objects.
 */
 /*
gsl_odeiv_step *
gsl_odeiv_step_new(
  const char * name,
  unsigned int dim,
  unsigned int ord,
  size_t state_size,
  size_t work_size);
*/


void gsl_odeiv_step_construct(
  gsl_odeiv_step *,
  const char * _name,
  int  (*_step)  (void *, double, double, double *, double *, const double *, double *, const gsl_odeiv_system *),
  int  (*_reset) (void *),
  void (*_free)  (void *),
  int can_use_dydt,
  size_t dimension,
  unsigned int order
  );


#endif /* !ODEIV_UTIL_H */
