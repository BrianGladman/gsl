/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


typedef  struct gsl_odeiv_step_rkck_struct  gsl_odeiv_step_rkck;

struct gsl_odeiv_step_rkck_struct
{
  gsl_odeiv_step parent;  /* inherits from gsl_odeiv_step */

  double * work;  /* generic work space */
};


static int  rkck_step(void * self, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
static void rkck_free(void * self);


gsl_odeiv_step *
gsl_odeiv_step_rkck_new(void)
{
  gsl_odeiv_step_rkck * s = (gsl_odeiv_step_rkck *) malloc(sizeof(gsl_odeiv_step_rkck));
  if(s != 0) {
    gsl_odeiv_step_construct(&(s->parent),
      "rkck",
      rkck_step,
      0,
      rkck_free,
      0,
      0,
      4);
    s->work = 0;
  }
  return (gsl_odeiv_step *) s;
}


static int
rkck_step(
  void * self,
  double t,
  double h,
  double y[],
  double yerr[],
  const double dydt_in[],
  double dydt_out[],
  const gsl_odeiv_system * sys)
{
  /* Cash-Karp constants */
  static const double ah[] = { 1.0/5.0, 0.3, 3.0/5.0, 1.0, 7.0/8.0 };
  static const double b21  = 1.0/5.0;
  static const double b3[] = { 3.0/40.0, 9.0/40.0 };
  static const double b4[] = { 0.3, -0.9, 1.2 };
  static const double b5[] = { -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0 };
  static const double b6[] = { 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0 };
  static const double c1 = 37.0/378.0;
  static const double c3 = 250.0/621.0;
  static const double c4 = 125.0/594.0;
  static const double c6 = 512.0/1771.0;
  const double ec[] = {
    0.0, c1 - 2825.0/27648.0, 0.0,
    c3 - 18575.0/48384.0, c4 - 13525.0/55296.0,
    -277.00/14336.0, c6 - 0.25
    };

  int i;
  int status = 0;
  size_t dim;

  double * k1;
  double * k2;
  double * k3;
  double * k4;
  double * k5;
  double * k6;
  double * ytmp;

  gsl_odeiv_step_rkck * my = (gsl_odeiv_step_rkck *) self;

  if(sys->dimension == 0) {
    return GSL_EINVAL;
  }

  if(sys->dimension != my->parent.dimension) {
    if(my->work != 0) free(my->work);
    my->parent.dimension = sys->dimension;
    my->work = (double *) malloc(7 * sys->dimension * sizeof(double));
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
  k4 = my->work + 3*dim;
  k5 = my->work + 4*dim;
  k6 = my->work + 5*dim;
  ytmp = my->work + 6*dim;

  /* k1 step */
  if(dydt_in != 0) {
    memcpy(k1, dydt_in, dim * sizeof(double));
  }
  else {
    status += ( GSL_ODEIV_FN_EVAL(sys, t, y, k1) != 0 );
  }
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + b21*h*k1[i];

  /* k2 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + ah[0]*h, ytmp, k2) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b3[0]*k1[i] + b3[1]*k2[i]);

  /* k3 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + ah[1]*h, ytmp, k3) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b4[0]*k1[i] + b4[1]*k2[i] + b4[2]*k3[i]);

  /* k4 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + ah[2]*h, ytmp, k4) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b5[0]*k1[i] + b5[1]*k2[i] + b5[2]*k3[i] + b5[3]*k4[i]);

  /* k5 step */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + ah[3]*h, ytmp, k5) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b6[0]*k1[i] + b6[1]*k2[i] + b6[2]*k3[i] + b6[3]*k4[i] + b6[4]*k5[i]);

  /* k6 step and final sum */
  status += ( GSL_ODEIV_FN_EVAL(sys, t + ah[4]*h, ytmp, k6) != 0 );
  for(i=0;i<dim;i++) {
    const double d_i = c1*k1[i] + c3*k3[i] + c4*k4[i] + c6*k6[i];
    y[i] += h * d_i;
    if(dydt_out != 0) dydt_out[i] = d_i;
  }

  /* difference between 4th and 5th order */
  for(i=0;i<dim;i++)
    yerr[i] = h*(ec[1]*k1[i] + ec[3]*k3[i] + ec[4]*k4[i] + ec[5]*k5[i] + ec[6]*k6[i]);

  return ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}


static void
rkck_free(void * self)
{
  if(self != 0) {
    gsl_odeiv_step_rkck * my = (gsl_odeiv_step_rkck *) self;
    if(my->work != 0) free(my->work);
    free(self);
  }
}
