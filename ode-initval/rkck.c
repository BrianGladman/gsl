/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include "odeiv_util.h"
#include "gsl_odeiv.h"


static gsl_odeiv_step * rkck_create(unsigned int dimension);
static int  rkck_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt);

const gsl_odeiv_step_factory gsl_odeiv_step_factory_rkck = 
{
  "rkck",
  rkck_create
};


static
gsl_odeiv_step *
rkck_create(unsigned int dimension)
{
  gsl_odeiv_step * step = gsl_odeiv_step_new(dimension, 4, 0, 7 * dimension * sizeof(double));
  step->_step = rkck_step;
  return step;
}


static
int
rkck_step(void * state, void * work, unsigned int dim, double t, double h, double y[], double yerr[], const gsl_odeiv_system * dydt)
{
  /* Cash-Karp constants */
  const double ah[] = { 1.0/5.0, 0.3, 3.0/5.0, 1.0, 7.0/8.0 };
  const double b21  = 1.0/5.0;
  const double b3[] = { 3.0/40.0, 9.0/40.0 };
  const double b4[] = { 0.3, -0.9, 1.2 };
  const double b5[] = { -11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0 };
  const double b6[] = { 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0 };
  const double c1 = 37.0/378.0;
  const double c3 = 250.0/621.0;
  const double c4 = 125.0/594.0;
  const double c6 = 512.0/1771.0;
  const double ec[] = {
    0.0, c1 - 2825.0/27648.0, 0.0,
    c3 - 18575.0/48384.0, c4 - 13525.0/55296.0,
    -277.00/14336.0, c6 - 0.25
    };

  int i;
  int status;

  /* divide up the work space */
  double * w  = (double *) work;
  double * k1 = w;
  double * k2 = w + dim;
  double * k3 = w + 2*dim;
  double * k4 = w + 3*dim;
  double * k5 = w + 4*dim;
  double * k6 = w + 5*dim;
  double * ytmp = w + 6*dim;

  /* k1 step */
  status  = ( GSL_ODEIV_FN_EVAL(dydt, t, y, k1) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + b21*h*k1[i];

  /* k2 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + ah[0]*h, ytmp, k2) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b3[0]*k1[i] + b3[1]*k2[i]);

  /* k3 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + ah[1]*h, ytmp, k3) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b4[0]*k1[i] + b4[1]*k2[i] + b4[2]*k3[i]);

  /* k4 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + ah[2]*h, ytmp, k4) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b5[0]*k1[i] + b5[1]*k2[i] + b5[2]*k3[i] + b5[3]*k4[i]);

  /* k5 step */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + ah[3]*h, ytmp, k5) != 0 );
  for(i=0;i<dim;i++)
    ytmp[i] = y[i] + h*(b6[0]*k1[i] + b6[1]*k2[i] + b6[2]*k3[i] + b6[3]*k4[i] + b6[4]*k5[i]);

  /* k6 step and final sum */
  status += ( GSL_ODEIV_FN_EVAL(dydt, t + ah[4]*h, ytmp, k6) != 0 );
  for(i=0;i<dim;i++)
    y[i] += h*(c1*k1[i] + c3*k3[i] + c4*k4[i] + c6*k6[i]);

  /* difference between 4th and 5th order */
  for(i=0;i<dim;i++)
    yerr[i] = h*(ec[1]*k1[i] + ec[3]*k3[i] + ec[4]*k4[i] + ec[5]*k5[i] + ec[6]*k6[i]);


  return ( status == 0 ? GSL_SUCCESS : GSL_EBADFUNC );
}
