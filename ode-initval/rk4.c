/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


/* rk4 state object */
typedef struct {
  double * work;   /* workspace  */
}
gsl_odeiv_step_rk4_state;


/* rk4 stepper object */
typedef struct {
  int  (*_step)  (void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
  int  (*_reset) (void * state);
  void (*_free)  (void * state);
  void * _state;
  unsigned int dimension;
}
gsl_odeiv_step_rk4;


static gsl_odeiv_step * rk4_create(unsigned int dimension);
static int  rk4_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
static int  rk4_reset(void *);
static void rk4_free(void *);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_rk4 = 
{
  "rk4",
  rk4_create
};


static
gsl_odeiv_step *
rk4_create(unsigned int dimension)
{
  gsl_odeiv_step_rk4 * rk4 = (gsl_odeiv_step_rk4 *) malloc(sizeof(gsl_odeiv_step_rk4));
  if(rk4_step != 0 && dimension > 0) {
    rk4->dimension = dimension;
    rk4->_step  = rk4_step;
    rk4->_free  = rk4_free;
    rk4->_reset = rk4_reset;
    rk4->_state = (gsl_odeiv_step_rk4_state *) malloc(sizeof(gsl_odeiv_step_rk4_state));
    if(rk4->_state != 0) {
      gsl_odeiv_step_rk4_state * s = (gsl_odeiv_step_rk4_state *) rk4->_state;
      s->work = (double *) malloc(3*dimension * sizeof(double));
      if(s->work == 0) {
        free(rk4->_state);
	free(rk4);
	return 0;
      }
    }
    else {
      free(rk4);
      return 0;
    }
  }
  return (gsl_odeiv_step *) rk4;
}


static
int
rk4_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt)
{
  int i;
  int status;
  gsl_odeiv_step_rk4_state * s = (gsl_odeiv_step_rk4_state *) state;

  /* divide up the workspace */
  double * k    = s->work;
  double * y0   = s->work + dim;
  double * ytmp = s->work + 2*dim;

  /* Copy the starting value. We will write over
   * the y[] vector, using it for scratch and
   * then filling it with the final result.
   */
  memcpy(y0, y, dim * sizeof(double));

  /* k1 step */
  status = ( GSL_ODEIV_FN_EVAL(dydt, t, y0, k) != 0 );
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
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + h, ytmp, k) != 0 );
  for(i=0; i<dim; i++) {
    y[i]   += h/6.0 * k[i];
    yerr[i] = h * y[i];
    y[i]   += y0[i];
  } 

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static
int
rk4_reset(void * state)
{
  /* gsl_odeiv_step_rk4_state * s = (gsl_odeiv_step_rk4_state *) state; */
  /* nothing to do for state reset */
  return GSL_SUCCESS;
}


static
void
rk4_free(void * state)
{
  gsl_odeiv_step_rk4_state * s = (gsl_odeiv_step_rk4_state *) state;
  if(s->work != 0) free(s->work);
  free(s);
}
