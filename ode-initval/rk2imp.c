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


typedef  struct gsl_odeiv_step_rk2imp_struct  gsl_odeiv_step_rk2imp;

struct gsl_odeiv_step_rk2imp_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  double * work;  /* generic work space */
};



static int  rk2imp_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static void rk2imp_free(void *);


gsl_odeiv_step *
gsl_odeiv_step_rk2imp_new(void)
{
  gsl_odeiv_step_rk2imp * s = (gsl_odeiv_step_rk2imp *) malloc(sizeof(gsl_odeiv_step_rk2imp));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "rk2imp",
      rk2imp_step,
      0,
      rk2imp_free,
      0,
      0,
      2);
    s->work = 0;
  }
  return (gsl_odeiv_step *) s;
}


static
int
rk2imp_step(
  void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  const int iter_steps = 3;
  int status = 0;
  int nu;
  int i;
  size_t dim;

  double * knu;
  double * ytmp;

  gsl_odeiv_step_rk2imp * my = (gsl_odeiv_step_rk2imp *) self;

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
  }

  dim = my->parent.dimension;

  /* divide up the workspace */
  knu  = my->work;
  ytmp = my->work + dim;

  /* initialization step */
  if(dydt_in != 0) {
    memcpy(knu, dydt_in, dim * sizeof(double));
  }
  else {
    status += ( GSL_ODEIV_FN_EVAL(sys, t, y, knu) != 0 );
  }

  /* iterative solution */
  for(nu=0; nu<iter_steps; nu++) {
    for(i=0; i<dim; i++) {
      ytmp[i] = y[i] + 0.5*h*knu[i];
    }
    status += ( GSL_ODEIV_FN_EVAL(sys, t + 0.5*h, ytmp, knu) != 0 );
  }

  /* assignment */
  for(i=0; i<dim; i++) {
    if(dydt_out != 0) dydt_out[i] = knu[i];
    y[i]   += h * knu[i];
    yerr[i] = h * h * knu[i];
  }

  return  ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static void
rk2imp_free(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_rk2imp * my = (gsl_odeiv_step_rk2imp *) self;
    if(my->work != 0) free(my->work);
    free(self);
  }
}
