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


/* gear2 state object */
typedef struct {
  int primed;                /* flag indicating that yim1 is ready */
  double last_h;             /* last step size */
  double * yim1;             /* y_{i-1}    */
  gsl_odeiv_step * primer;   /* stepper to use for priming */
}
gear2_state;


static gsl_odeiv_step * gear2_create(unsigned int dimension);
static int  gear2_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static int  gear2_reset(void * state);
static void gear2_free(void * state, void * work);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_gear2 = 
{
  "gear2",
  gear2_create
};


static
gsl_odeiv_step *
gear2_create(unsigned int dimension)
{
  gear2_state * state;

  gsl_odeiv_step * step = gsl_odeiv_step_new(gsl_odeiv_step_factory_gear2.name, dimension, 3, sizeof(gear2_state), 2*dimension * sizeof(double));
  if(step == 0) return 0;

  step->_step  = gear2_step;
  step->_reset = gear2_reset;
  step->_free  = gear2_free;
  step->can_use_dydt = 0;
  step->stutter = 0;

  state = (gear2_state *) step->_state;
  state->primed = 0;
  state->last_h = 0.0;
  state->yim1 = (double *)malloc(dimension * sizeof(double));
  state->primer = gsl_odeiv_step_factory_rk4imp.create(dimension);

  if(state->primer == 0 || state->yim1 == 0) {
    gsl_odeiv_step_free(step);
    return 0;
  }
  
  return step;
}


static
int
gear2_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt)
{
  gsl_odeiv_step * my = (gsl_odeiv_step *) self;

  gear2_state * s = (gear2_state *) my->_state;

  const unsigned int dim = my->dimension;

  my->stutter = 0;

  if(s->primed && fabs((s->last_h - h)/h) < 8.0 * GSL_DBL_EPSILON) {
    /* We have a previous y value in the buffer, and the step
     * sizes match, so we go ahead with the Gear step.
     */

    /* divide up the workspace */
    double * w  = (double *) my->_work;
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
        y[i] = ((4.0 * y0[i] - s->yim1[i]) + 2.0 * h * k[i]) / 3.0;
      }
    }

    /* Estimate error and update the state buffer. */
    for(i=0; i<dim; i++) {
      yerr[i]    = h * h * (y[i] - y0[i]); /* error is third order */
      s->yim1[i] = y0[i];
    }

    /* Make note of step size. */
    s->last_h = h;

    return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
  }
  else {
    /* Execute a single-step method to prime the process.
     * Note that we do this if the step size changes, so
     * frequent step size changes will cause the method
     * to stutter.
     */
    int status;
    memcpy(s->yim1, y, dim * sizeof(double));
    status = gsl_odeiv_step_impl(s->primer, t, h, y, yerr, dydt_in, dydt_out, dydt);

    /* Make note of step size and indicate readiness for a Gear step. */
    s->last_h = h;
    s->primed = 1;
    my->stutter = 1;

    return status;
  }
}


static
int
gear2_reset(void * state)
{
  gear2_state * s = (gear2_state *) state;
  s->primed = 0;
  s->last_h = 0.0;
  return GSL_SUCCESS;
}


static
void
gear2_free(void * state, void * work)
{
  gear2_state * s = (gear2_state *) state;
  if(s != 0) {
    if(s->yim1 != 0) free(s->yim1);
    gsl_odeiv_step_free(s->primer);
  }
}
