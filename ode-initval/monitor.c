/* ode-initval/monitor.c
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
#include "gsl_odeiv.h"


static int
dump_step(void * self, double t, unsigned int dim, const double y[], const double yerr[])
{
  unsigned int i;
  gsl_odeiv_evolve_mon * m = (gsl_odeiv_evolve_mon *) self;
  fprintf(m->f, "%20.16g", t);
  for(i=0; i<dim; i++) {
    fprintf(m->f, "  %22.18g %22.18g", y[i], yerr[i]);
  }
  fprintf(m->f, "\n");
  return GSL_SUCCESS;
}


gsl_odeiv_evolve_mon *
gsl_odeiv_evolve_mon_stream_new(FILE * f_in)
{
  gsl_odeiv_evolve_mon * m = (gsl_odeiv_evolve_mon *) malloc(sizeof(gsl_odeiv_evolve_mon));
  if(m != 0) {
    m->f = f_in;
    m->pre_step = 0;
    m->post_step = dump_step;
    m->params = 0;
  }
  return m;
}


void
gsl_odeiv_evolve_mon_free(gsl_odeiv_evolve_mon * mon)
{
  if(mon != 0) {
    if(mon->params != 0) free(mon->params);
    free(mon);
  }
}
