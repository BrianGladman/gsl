/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


/* rk2 state object */
typedef struct {
  double * work;   /* workspace  */
}
gsl_odeiv_step_rk2_state;


/* rk2 stepper object */
/*
typedef struct {
  int  (*_step)  (void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_system * dydt);
  int  (*_reset) (void * state);
  void (*_free)  (void * state);
  void * _state;
  unsigned int dimension;
}
gsl_odeiv_step_rk2;
*/
typedef gsl_odeiv_step gsl_odeiv_step_rk2;


static gsl_odeiv_step * rk2_create(unsigned int dimension);
static int  rk2_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt);
static int  rk2_reset(void *);
static void rk2_free(void *);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_rk2 = 
{
  "rk2",
  rk2_create
};


static
gsl_odeiv_step *
rk2_create(unsigned int dimension)
{
  gsl_odeiv_step_rk2 * rk2;

  if(dimension == 0) return 0;

  rk2 = (gsl_odeiv_step_rk2 *) malloc(sizeof(gsl_odeiv_step_rk2));
  if(rk2_step != 0) {
    rk2->dimension = dimension;
    rk2->_step  = rk2_step;
    rk2->_free  = rk2_free;
    rk2->_reset = rk2_reset;
    rk2->_state = (gsl_odeiv_step_rk2_state *) malloc(sizeof(gsl_odeiv_step_rk2_state));
    if(rk2->_state != 0) {
      gsl_odeiv_step_rk2_state * s = (gsl_odeiv_step_rk2_state *) rk2->_state;
      s->work = (double *) malloc(4*dimension * sizeof(double));
      if(s->work == 0) {
        free(rk2->_state);
	free(rk2);
	return 0;
      }
    }
    else {
      free(rk2);
      return 0;
    }
  }
  return (gsl_odeiv_step *) rk2;
}


static
int
rk2_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt)
{
  gsl_odeiv_step_rk2_state * s = (gsl_odeiv_step_rk2_state *) state;

  int i;
  int status;

  /* divide up the work space */
  double * k1 = s->work;
  double * k2 = s->work + dim;
  double * k3 = s->work + 2*dim;
  double * ytmp = s->work + 3*dim;

  /* k1 step */
  status  = ( GSL_ODEIV_FN_EVAL(dydt, t, y, k1) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + 0.5*h*k1[i];

  /* k2 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h, ytmp, k2) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(-k1[i] + 2.0*k2[i]);

  /* k3 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + h, ytmp, k3) != 0 );

  /* final sum and error estimate */
  for(i=0;i<dim;i++) {
    const double ksum3 = (k1[i] + 4.0*k2[i] + k3[i])/6.0;
    y[i]    += h*ksum3;
    yerr[i]  = h*(k2[i] - ksum3);
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static
int
rk2_reset(void * state)
{
  /* gsl_odeiv_step_rk2_state * s = (gsl_odeiv_step_rk2_state *) state; */
  /* nothing to do for state reset */
  return GSL_SUCCESS;
}


static
void
rk2_free(void * state)
{
  gsl_odeiv_step_rk2_state * s = (gsl_odeiv_step_rk2_state *) state;
  if(s->work != 0) free(s->work);
  free(s);
}
