/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


/* gear2 state object */
typedef struct {
  double * work;             /* workspace  */
  double * yim1;             /* y_{i-1}    */
  int primed;                /* flag indicating that yim1 is ready */
  gsl_odeiv_step * primer;   /* stepper to use for priming */
  double last_h;             /* last step size */
}
gsl_odeiv_step_gear2_state;


/* gear2 stepper object */
typedef struct {
  int  (*_step)  (void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
  int  (*_reset) (void * state);
  void (*_free)  (void * state);
  void * _state;
  unsigned int dimension;
}
gsl_odeiv_step_gear2;


static gsl_odeiv_step * gear2_create(unsigned int dimension);
static int  gear2_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
static int  gear2_reset(void *);
static void gear2_free(void *);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_gear2 = 
{
  "gear2",
  gear2_create
};


static
gsl_odeiv_step *
gear2_create(unsigned int dimension)
{
  gsl_odeiv_step_gear2 * gear2;

  if(dimension == 0) return 0;

  gear2 = (gsl_odeiv_step_gear2 *) malloc(sizeof(gsl_odeiv_step_gear2));
  if(gear2_step != 0) {
    gear2->dimension = dimension;
    gear2->_step  = gear2_step;
    gear2->_free  = gear2_free;
    gear2->_reset = gear2_reset;
    gear2->_state = (gsl_odeiv_step_gear2_state *) malloc(sizeof(gsl_odeiv_step_gear2_state));
    if(gear2->_state != 0) {
      gsl_odeiv_step_gear2_state * s = (gsl_odeiv_step_gear2_state *) gear2->_state;
      s->work = (double *) malloc(2*dimension * sizeof(double));
      if(s->work == 0) {
        free(gear2->_state);
	free(gear2);
	return 0;
      }
      s->primer = gsl_odeiv_step_factory_rk4imp.create(dimension);
      if(s->primer == 0) {
        free(s->work);
        free(gear2->_state);
	free(gear2);
      }
      s->yim1 = (double *) malloc(dimension * sizeof(double));
      if(s->yim1 == 0) {
        gsl_odeiv_step_free(s->primer);
        free(s->work);
        free(gear2->_state);
	free(gear2);
	return 0;
      }
    }
    else {
      free(gear2);
      return 0;
    }
  }
  return (gsl_odeiv_step *) gear2;
}


static
int
gear2_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt)
{
  gsl_odeiv_step_gear2_state * s = (gsl_odeiv_step_gear2_state *) state;

  if(s->primed && fabs((s->last_h - h)/h) < 8.0 * GSL_DBL_EPSILON) {
    /* We have a previous y value in the buffer, and the step
     * sizes match, so we go ahead with the Gear step.
     */

    /* divide up the workspace */
    double * k = s->work;
    double * y0 = s->work + dim;

    const int iter_steps = 3;
    int status = 0;
    int nu;
    int i;

    memcpy(y0, y, dim * sizeof(double));

    /* iterative solution */
    for(nu=0; nu<iter_steps; nu++) {
      status += ( GSL_ODEIV_FN_EVAL(dydt, t + h, y, k) != 0 );
      for(i=0; i<dim; i++) {
        y[i] = ((4.0 * y0[i] - s->yim1[i]) + 2.0 * h * k[i]) / 3.0;
      }
    }

    /* Estimate error and update the buffer. */
    for(i=0; i<dim; i++) {
      yerr[i] = h * (y[i] - y0[i]);
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
    int status = gsl_odeiv_step_impl(s->primer, t, h, y, yerr, dydt);
    s->primed = 1;

    /* Make note of step size. */
    s->last_h = h;

    return status;
  }
}


static
int
gear2_reset(void * state)
{
  gsl_odeiv_step_gear2_state * s = (gsl_odeiv_step_gear2_state *) state;
  s->primed = 0;
  s->last_h = 0.0;
  return GSL_SUCCESS;
}


static
void
gear2_free(void * state)
{
  gsl_odeiv_step_gear2_state * s = (gsl_odeiv_step_gear2_state *) state;
  if(s->work != 0) free(s->work);
  if(s->yim1 != 0) free(s->yim1);
  gsl_odeiv_step_free(s->primer);
  free(s);
}
