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


static gsl_odeiv_step * rk4imp_create(unsigned int dimension);
static int rk4imp_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_rk4imp = 
{
  "rk4imp",
  rk4imp_create
};


static
gsl_odeiv_step *
rk4imp_create(unsigned int dimension)
{
  gsl_odeiv_step * step = gsl_odeiv_step_new(gsl_odeiv_step_factory_rk4imp.name, dimension, 4, 0, 4*dimension * sizeof(double));
  step->_step = rk4imp_step;
  return step;
}


static
int
rk4imp_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt)
{
  /* divide up the workspace */
  double * w    = (double *) work;
  double * k1nu = w;
  double * k2nu = w + dim;
  double * ytmp1 = w + 2*dim;
  double * ytmp2 = w + 3*dim;

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
