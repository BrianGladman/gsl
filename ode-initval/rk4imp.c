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


typedef  struct gsl_odeiv_step_rk4imp_struct  gsl_odeiv_step_rk4imp;

struct gsl_odeiv_step_rk4imp_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  double * work;  /* generic work space */
};


static int  rk4imp_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static void rk4imp_free(void * self);


gsl_odeiv_step *
gsl_odeiv_step_rk4imp_new(void)
{
  gsl_odeiv_step_rk4imp * s = (gsl_odeiv_step_rk4imp *) malloc(sizeof(gsl_odeiv_step_rk4imp));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "rk4imp",
      rk4imp_step,
      0,
      rk4imp_free,
      0,
      0,
      2);
    s->work = 0;
  }
  return (gsl_odeiv_step *) s;
}


static int
rk4imp_step(
  void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  const double ir3 = 1.0/M_SQRT3;
  const int iter_steps = 3;
  int status = 0;
  int nu;
  int i;
  size_t dim;
  
  double * k1nu;
  double * k2nu;
  double * ytmp1;
  double * ytmp2;

  gsl_odeiv_step_rk4imp * my = (gsl_odeiv_step_rk4imp *) self;

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

  /* divide up the workspace */
  k1nu = my->work;
  k2nu = my->work + dim;
  ytmp1 = my->work + 2*dim;
  ytmp2 = my->work + 3*dim;

  /* initialization step */
  if(dydt_in != 0) {
    memcpy(k1nu, dydt_in, dim * sizeof(double));
  }
  else {
    status += ( GSL_ODEIV_FN_EVAL(sys, t, y, k1nu) != 0 );
  }
  memcpy(k2nu, k1nu, dim * sizeof(double));

  /* iterative solution */
  for(nu=0; nu<iter_steps; nu++) {
    for(i=0; i<dim; i++) {
      ytmp1[i] = y[i] + h*(0.25 * k1nu[i] + 0.5*(0.5 - ir3) * k2nu[i]);
      ytmp2[i] = y[i] + h*(0.25 * k2nu[i] + 0.5*(0.5 + ir3) * k1nu[i]);
    }
    status += ( GSL_ODEIV_FN_EVAL(sys, t + 0.5*h*(1.0-ir3), ytmp1, k1nu) != 0 );
    status += ( GSL_ODEIV_FN_EVAL(sys, t + 0.5*h*(1.0+ir3), ytmp2, k2nu) != 0 );
  }

  /* assignment */
  for(i=0; i<dim; i++) {
    const double d_i = 0.5 * (k1nu[i] + k2nu[i]);
    if(dydt_out != 0) dydt_out[i] = d_i;
    y[i]   += h * d_i;
    yerr[i] = h * h * d_i; /* FIXME: is this an overestimate ? */
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static void
rk4imp_free(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_rk4imp * my = (gsl_odeiv_step_rk4imp *) self;
    if(my->work != 0) free(my->work);
    free(self);
  }
}
