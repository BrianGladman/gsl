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
  double * y0;    /* last value */
  double * k;     /* workspace  */
  double * ytmp;  /* workspace */
}
gsl_odeiv_step_rk4_state;


/* rk4 stepper object */
typedef struct {
  int  (*_step)  (void * state, unsigned int dim, double t, double h, double y[], gsl_odeiv_function * dydt);
  int  (*_reset) (void * state);
  void (*_free)  (void * state);
  void * _state;
  unsigned int dimension;
}
gsl_odeiv_step_rk4;


static gsl_odeiv_step * rk4_create(unsigned int dimension);
static int  rk4_step(void * state, unsigned int dim, double t, double h, double y[], gsl_odeiv_function * dydt);
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
      s->k = (double *) malloc(dimension * sizeof(double));
      if(s->k == 0) {
        free(rk4->_state);
	free(rk4);
	return 0;
      }
      s->y0 = (double *) malloc(dimension * sizeof(double));
      if(s->y0 == 0) {
        free(s->k);
        free(rk4->_state);
	free(rk4);
	return 0;
      }
      s->ytmp = (double *) malloc(dimension * sizeof(double));
      if(s->ytmp == 0) {
        free(s->y0);
        free(s->k);
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
rk4_step(void * state, unsigned int dim, double t, double h, double y[], gsl_odeiv_function * dydt)
{
  int i;
  int status;
  gsl_odeiv_step_rk4_state * s = (gsl_odeiv_step_rk4_state *) state;

  memcpy(s->y0, y, dim * sizeof(double));

  /* k1 step */
  status = ( GSL_ODEIV_FN_EVAL(dydt, t, s->y0, s->k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/6.0 * s->k[i];
    s->ytmp[i] = s->y0[i] + 0.5*h * s->k[i];
  }

  /* k2 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h, s->ytmp, s->k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/3.0 * s->k[i];
    s->ytmp[i] = s->y0[i] + 0.5*h * s->k[i];
  }  

  /* k3 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h, s->ytmp, s->k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/3.0 * s->k[i];
    s->ytmp[i] = s->y0[i] + h * s->k[i];
  } 

  /* k4 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + h, s->ytmp, s->k) != 0 );
  for(i=0; i<dim; i++) {
    y[i] += h/6.0 * s->k[i];
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
  if(s->k != 0) free(s->k);
  free(s);
}
