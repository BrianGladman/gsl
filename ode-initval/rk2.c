/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


static gsl_odeiv_step * rk2_create(unsigned int dimension);
static int  rk2_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt);


const gsl_odeiv_step_factory gsl_odeiv_step_factory_rk2 = 
{
  "rk2",
  rk2_create
};


static
gsl_odeiv_step *
rk2_create(unsigned int dimension)
{
  gsl_odeiv_step * step = gsl_odeiv_step_new(dimension, 2, 0, 4*dimension * sizeof(double));
  step->_step = rk2_step;
  return step;
}


static
int
rk2_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt)
{
  int i;
  int status;

  /* divide up the work space */
  double * w  = (double *) work;
  double * k1 = w;
  double * k2 = w + dim;
  double * k3 = w + 2*dim;
  double * ytmp = w + 3*dim;

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
