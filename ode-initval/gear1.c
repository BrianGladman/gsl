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


static gsl_odeiv_step * gear1_create(unsigned int dimension);
static int gear1_step(void *, void *, unsigned int dim, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_gear1 = 
{
  "gear1",
  gear1_create
};


static
gsl_odeiv_step *
gear1_create(unsigned int dimension)
{
  gsl_odeiv_step * step = gsl_odeiv_step_new(gsl_odeiv_step_factory_gear1.name, dimension, 1 /* FIXME: ?? */, 0, 2*dimension * sizeof(double));
  step->_step = gear1_step;
  step->can_use_dydt = 0;
  return step;
}


static
int
gear1_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt)
{
  /* divide up the workspace */
  double * w  = (double *) work;
  double * k  = w;
  double * y0 = w + dim;

  const int iter_steps = 3;
  int status = 0;
  int nu;
  int i;

  memcpy(y0, y, dim * sizeof(double));

  /* iterative solution */
  if(dydt_out != 0) {
    k = dydt_out;
  }
  for(nu=0; nu<iter_steps; nu++) {
    status += ( GSL_ODEIV_FN_EVAL(dydt, t + h, y, k) != 0 );
    for(i=0; i<dim; i++) {
      y[i] = y0[i] + h * k[i];
    }
  }

  /* fudge the error estimate */
  for(i=0; i<dim; i++) {
    yerr[i] = h * h * k[i];
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}

