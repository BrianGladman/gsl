/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


typedef  struct gsl_odeiv_step_rk2_struct  gsl_odeiv_step_rk2;

struct gsl_odeiv_step_rk2_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  double * work;  /* generic work space */
};


static int  rk2_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static void rk2_free(void * self);



gsl_odeiv_step *
gsl_odeiv_step_rk2_new(void)
{
  gsl_odeiv_step_rk2 * s = (gsl_odeiv_step_rk2 *) malloc(sizeof(gsl_odeiv_step_rk2));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "rk2",
      rk2_step,
      0,
      rk2_free,
      0,
      0,
      2);
    s->work = 0;
  }
  return (gsl_odeiv_step *) s;
}


static int
rk2_step(
  void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  int i;
  int status = 0;
  size_t dim;
  double * k1;
  double * k2;
  double * k3;
  double * ytmp;

  gsl_odeiv_step_rk2 * my = (gsl_odeiv_step_rk2 *) self;

  if(sys->dimension == 0) {
    return GSL_EINVAL;
  }

  if(sys->dimension != my->parent.dimension) {
    if(my->work != 0) free(my->work);
    my->parent.dimension = sys->dimension;
    my->work = (double *) malloc(4 * sys->dimension * sizeof(double));
    if(my->work == 0) {
      my->parent.dimension = 0;
      return GSL_ENOMEM;
    }
  }

  dim = my->parent.dimension;

  /* divide up the work space */
  k1 = my->work;
  k2 = my->work + dim;
  k3 = my->work + 2*dim;
  ytmp = my->work + 3*dim;

  /* k1 step */
  if(dydt_in != 0) {
    memcpy(k1, dydt_in, dim * sizeof(double));
  }
  else {
    status += ( GSL_ODEIV_FN_EVAL(sys, t, y, k1) != 0 );
  }
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + 0.5*h*k1[i];

  /* k2 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + 0.5*h, ytmp, k2) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(-k1[i] + 2.0*k2[i]);

  /* k3 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + h, ytmp, k3) != 0 );

  /* final sum and error estimate */
  for(i=0;i<dim;i++) {
    const double ksum3 = (k1[i] + 4.0*k2[i] + k3[i])/6.0;
    if(dydt_out != 0) dydt_out[i] = ksum3;
    y[i]    += h*ksum3;
    yerr[i]  = h*(k2[i] - ksum3);
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static void
rk2_free(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_rk2 * my = (gsl_odeiv_step_rk2 *) self;
    if(my->work != 0) free(my->work);
    free(self);
  }
}
