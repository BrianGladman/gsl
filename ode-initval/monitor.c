/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl/gsl_errno.h>
#include "gsl_odeiv.h"


static int
dump_step(void * self, double t, unsigned int dim, const double y[], const double yerr[])
{
  int i;
  gsl_odeiv_evolve_mon * m = (gsl_odeiv_evolve_mon *) self;
  fprintf(m->f, "%20.16g", t);
  for(i=0; i<dim; i++) {
    fprintf(m->f, "  %22.18g", y[i]);
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
