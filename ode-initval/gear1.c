/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


/* gear1 state object */
typedef struct {
  double * work;   /* workspace  */
}
gsl_odeiv_step_gear1_state;


/* gear1 stepper object */
typedef struct {
  int  (*_step)  (void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_system * dydt);
  int  (*_reset) (void * state);
  void (*_free)  (void * state);
  void * _state;
  unsigned int dimension;
}
gsl_odeiv_step_gear1;


static gsl_odeiv_step * gear1_create(unsigned int dimension);
static int  gear1_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_system * dydt);
static int  gear1_reset(void *);
static void gear1_free(void *);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_gear1 = 
{
  "gear1",
  gear1_create
};


static
gsl_odeiv_step *
gear1_create(unsigned int dimension)
{
  gsl_odeiv_step_gear1 * gear1;

  if(dimension == 0) return 0;

  gear1 = (gsl_odeiv_step_gear1 *) malloc(sizeof(gsl_odeiv_step_gear1));
  if(gear1_step != 0) {
    gear1->dimension = dimension;
    gear1->_step  = gear1_step;
    gear1->_free  = gear1_free;
    gear1->_reset = gear1_reset;
    gear1->_state = (gsl_odeiv_step_gear1_state *) malloc(sizeof(gsl_odeiv_step_gear1_state));
    if(gear1->_state != 0) {
      gsl_odeiv_step_gear1_state * s = (gsl_odeiv_step_gear1_state *) gear1->_state;
      s->work = (double *) malloc(2*dimension * sizeof(double));
      if(s->work == 0) {
        free(gear1->_state);
	free(gear1);
	return 0;
      }
    }
    else {
      free(gear1);
      return 0;
    }
  }
  return (gsl_odeiv_step *) gear1;
}


static
int
gear1_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_system * dydt)
{
  gsl_odeiv_step_gear1_state * s = (gsl_odeiv_step_gear1_state *) state;

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
      y[i] = y0[i] + h * k[i];
    }
  }

  for(i=0; i<dim; i++) {
    yerr[i] = h * h * k[i];
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static
int
gear1_reset(void * state)
{
  /* nothing to do for state reset */
  return GSL_SUCCESS;
}


static
void
gear1_free(void * state)
{
  gsl_odeiv_step_gear1_state * s = (gsl_odeiv_step_gear1_state *) state;
  if(s->work != 0) free(s->work);
  free(s);
}
