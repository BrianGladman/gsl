/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


static gsl_odeiv_step * rk4_create(unsigned int dimension);
static int rk4_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_rk4 = 
{
  "rk4",
  rk4_create
};


static
gsl_odeiv_step *
rk4_create(unsigned int dimension)
{
  gsl_odeiv_step * step = gsl_odeiv_step_new(gsl_odeiv_step_factory_rk4.name, dimension, 4, 0, 3*dimension * sizeof(double));
  step->_step = rk4_step;
  step->can_use_dydt = 1;
  return step;
}


static
int
rk4_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt)
{
  /* divide up the workspace */
  double * w    = (double *) work;
  double * k    = w;
  double * y0   = w + dim;
  double * ytmp = w + 2*dim;

  int i;
  int status = 0;

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
    status += ( GSL_ODEIV_FN_EVAL(dydt, t, y0, k) != 0 );
  }
  for(i=0; i<dim; i++) {
    y[i] = h/6.0 * k[i]; /* use y[] to store delta_y */
    ytmp[i] = y0[i] + 0.5*h * k[i];
  }

  /* k2 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h, ytmp, k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/3.0 * k[i];
    ytmp[i] = y0[i] + 0.5*h * k[i];
  }  

  /* k3 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h, ytmp, k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/3.0 * k[i];
    ytmp[i] = y0[i] + h * k[i];
  } 

  /* k4 step, error estimate, and final sum */
  if(dydt_out != 0) {
    k = dydt_out;
  }
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + h, ytmp, k) != 0 );
  for(i=0; i<dim; i++) {
    y[i]   += h/6.0 * k[i];
    yerr[i] = h * y[i];
    y[i]   += y0[i];
  } 

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}
