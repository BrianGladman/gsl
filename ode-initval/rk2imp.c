/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


static gsl_odeiv_step * rk2imp_create(unsigned int dimension);
static int rk2imp_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_rk2imp = 
{
  "rk2imp",
  rk2imp_create
};


static
gsl_odeiv_step *
rk2imp_create(unsigned int dimension)
{
  gsl_odeiv_step * step = gsl_odeiv_step_new(gsl_odeiv_step_factory_rk2imp.name, dimension, 2, 0, 2*dimension * sizeof(double));
  step->_step = rk2imp_step;
  step->can_use_dydt = 1;
  step->stutter = 0;
  return step;
}


static
int
rk2imp_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt)
{
  gsl_odeiv_step * my = (gsl_odeiv_step *) self;

  const unsigned int dim = my->dimension;

  /* divide up the workspace */
  double * w    = (double *) my->_work;
  double * knu  = w;
  double * ytmp = w + dim;

  const int iter_steps = 3;
  int status = 0;
  int nu;
  int i;

  /* initialization step */
  if(dydt_in != 0) {
    memcpy(knu, dydt_in, dim * sizeof(double));
  }
  else {
    status += ( GSL_ODEIV_FN_EVAL(dydt, t, y, knu) != 0 );
  }

  /* iterative solution */
  for(nu=0; nu<iter_steps; nu++) {
    for(i=0; i<dim; i++) {
      ytmp[i] = y[i] + 0.5*h*knu[i];
    }
    status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h, ytmp, knu) != 0 );
  }

  /* assignment */
  for(i=0; i<dim; i++) {
    if(dydt_out != 0) dydt_out[i] = knu[i];
    y[i]   += h * knu[i];
    yerr[i] = h * h * knu[i];
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}
