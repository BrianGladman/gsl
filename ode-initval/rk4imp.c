/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_odeiv.h"


/* rk4imp state object */
typedef struct {
  double * work;   /* workspace  */
}
gsl_odeiv_step_rk4imp_state;


/* rk4imp stepper object */
typedef struct {
  int  (*_step)  (void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
  int  (*_reset) (void * state);
  void (*_free)  (void * state);
  void * _state;
  unsigned int dimension;
}
gsl_odeiv_step_rk4imp;


static gsl_odeiv_step * rk4imp_create(unsigned int dimension);
static int  rk4imp_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt);
static int  rk4imp_reset(void *);
static void rk4imp_free(void *);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_rk4imp = 
{
  "rk4imp",
  rk4imp_create
};


static
gsl_odeiv_step *
rk4imp_create(unsigned int dimension)
{
  gsl_odeiv_step_rk4imp * rk4imp = (gsl_odeiv_step_rk4imp *) malloc(sizeof(gsl_odeiv_step_rk4imp));
  if(rk4imp_step != 0 && dimension > 0) {
    rk4imp->dimension = dimension;
    rk4imp->_step  = rk4imp_step;
    rk4imp->_free  = rk4imp_free;
    rk4imp->_reset = rk4imp_reset;
    rk4imp->_state = (gsl_odeiv_step_rk4imp_state *) malloc(sizeof(gsl_odeiv_step_rk4imp_state));
    if(rk4imp->_state != 0) {
      gsl_odeiv_step_rk4imp_state * s = (gsl_odeiv_step_rk4imp_state *) rk4imp->_state;
      s->work = (double *) malloc(4*dimension * sizeof(double));
      if(s->work == 0) {
        free(rk4imp->_state);
	free(rk4imp);
	return 0;
      }
    }
    else {
      free(rk4imp);
      return 0;
    }
  }
  return (gsl_odeiv_step *) rk4imp;
}


static
int
rk4imp_step(void * state, unsigned int dim, double t, double h, double y[], double yerr[], gsl_odeiv_function * dydt)
{
  gsl_odeiv_step_rk4imp_state * s = (gsl_odeiv_step_rk4imp_state *) state;

  /* divide up the workspace */
  double * k1nu = s->work;
  double * k2nu = s->work + dim;
  double * ytmp1 = s->work + 2*dim;
  double * ytmp2 = s->work + 3*dim;

  const double ir3 = 1.0/M_SQRT3;
  const int iter_steps = 3;
  int status = 0;
  int nu;
  int i;

  /* initialization step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t, y, k1nu) != 0 );
  memcpy(k2nu, k1nu, dim * sizeof(double));

  /* iterative solution */
  for(nu=0; nu<iter_steps; nu++) {
    for(i=0; i<dim; i++) {
      ytmp1[i] = y[i] + h*(0.25 * k1nu[i] + 0.5*(0.5 - ir3) * k2nu[i]);
      ytmp2[i] = y[i] + h*(0.25 * k2nu[i] + 0.5*(0.5 + ir3) * k1nu[i]);
    }
    status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h*(1.0-ir3), ytmp1, k1nu) != 0 );
    status += ( GSL_ODEIV_FN_EVAL(dydt, t + 0.5*h*(1.0+ir3), ytmp2, k2nu) != 0 );
  }

  /* assignment */
  for(i=0; i<dim; i++) {
    const double del_i = 0.5 * h * (k1nu[i] + k2nu[i]);
    y[i]   += del_i;
    yerr[i] = h * del_i;
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static
int
rk4imp_reset(void * state)
{
  /* nothing to do for state reset */
  return GSL_SUCCESS;
}


static
void
rk4imp_free(void * state)
{
  gsl_odeiv_step_rk4imp_state * s = (gsl_odeiv_step_rk4imp_state *) state;
  if(s->work != 0) free(s->work);
  free(s);
}
