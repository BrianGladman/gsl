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


typedef  struct gsl_odeiv_step_gear2_struct  gsl_odeiv_step_gear2;

struct gsl_odeiv_step_gear2_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */
  gear2_state state;      /* state information  */

  double * work;          /* generic work space */

  int stutter;
};


static int  gear2_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static int  gear2_reset(void * self);
static void gear2_free(void * self);


gsl_odeiv_step *
gsl_odeiv_step_gear2_new(void)
{
  gsl_odeiv_step_gear2 * s = (gsl_odeiv_step_gear2 *) malloc(sizeof(gsl_odeiv_step_gear2));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "gear2",
      gear2_step,
      gear2_reset,
      gear2_free,
      0,
      0,
      3);
    s->work = 0;
    s->stutter = 0;
    s->state.yim1 = 0;
    s->state.primer = 0;
  }
  return (gsl_odeiv_step *) s;
}


static int
gear2_step(
  void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  const unsigned int dim = sys->dimension;

  gsl_odeiv_step_gear2 * my = (gsl_odeiv_step_gear2 *) self;

  if(sys->dimension == 0) {
    return GSL_EINVAL;
  }

  if(sys->dimension != my->parent.dimension) {

    if(my->work != 0) free(my->work);
    my->parent.dimension = sys->dimension;
    my->work = (double *) malloc(2 * sys->dimension * sizeof(double));
    if(my->work == 0) {
      my->parent.dimension = 0;
      return GSL_ENOMEM;
    }

    if(my->state.yim1 != 0) free(my->state.yim1);
    my->state.yim1 = (double *) malloc(sys->dimension * sizeof(double));
    if(my->state.yim1 == 0) {
      free(my->work);
      my->parent.dimension = 0;
      return GSL_ENOMEM;
    }

    if(my->state.primer == 0) {
      my->state.primer = gsl_odeiv_step_rk4imp_new();
      if(my->state.primer == 0) {
        free(my->work);
        free(my->state.yim1);
        return GSL_ENOMEM;
      }
    }

    my->state.primed = 0;
    my->state.last_h = 0.0;
  }

  my->stutter = 0;

  if(my->state.primed && fabs((my->state.last_h - h)/h) < 8.0 * GSL_DBL_EPSILON) {
    /* We have a previous y value in the buffer, and the step
     * sizes match, so we go ahead with the Gear step.
     */

    double * k  = my->work;
    double * y0 = my->work + dim;

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
      status += ( GSL_ODEIV_FN_EVAL(sys, t + h, y, k) != 0 );
      for(i=0; i<dim; i++) {
        y[i] = ((4.0 * y0[i] - my->state.yim1[i]) + 2.0 * h * k[i]) / 3.0;
      }
    }

    /* Estimate error and update the state buffer. */
    for(i=0; i<dim; i++) {
      yerr[i]    = h * h * (y[i] - y0[i]); /* error is third order */
      my->state.yim1[i] = y0[i];
    }

    /* Make note of step size. */
    my->state.last_h = h;

    return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
  }
  else {
    /* Execute a single-step method to prime the process.
     * Note that we do this if the step size changes, so
     * frequent step size changes will cause the method
     * to stutter.
     */
    int status;
    memcpy(my->state.yim1, y, dim * sizeof(double));
    status = gsl_odeiv_step_impl(my->state.primer, t, h, y, yerr, dydt_in, dydt_out, sys);

    /* Make note of step size and indicate readiness for a Gear step. */
    my->state.last_h = h;
    my->state.primed = 1;
    my->stutter = 1;

    return status;
  }
}


static int
gear2_reset(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_gear2 * my = (gsl_odeiv_step_gear2 *) self;
    my->state.primed = 0;
    my->state.last_h = 0.0;
    return GSL_SUCCESS;
  }
  else {
    return GSL_EFAULT;
  }
}


static void
gear2_free(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_gear2 * my = (gsl_odeiv_step_gear2 *) self;
    if(my->work != 0) free(my->work);
    if(my->state.yim1 != 0) free(my->state.yim1);
    gsl_odeiv_step_free(my->state.primer);
  }
}
