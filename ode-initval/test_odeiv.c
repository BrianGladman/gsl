/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl_test.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_ieee_utils.h>
#include "gsl_odeiv.h"


/* RHS for a + b t */

int rhs_linear(double t, const double y[], double f[], void * params)
{
  f[0] = 0.0;
  f[1] = y[0];
  return GSL_SUCCESS;
}

gsl_odeiv_function rhs_func_lin = {
  rhs_linear,
  0
};


/* RHS for sin(t),cos(t) */

int rhs_sin(double t, const double y[], double f[], void * params)
{
  f[0] = -y[1];
  f[1] =  y[0];
  return GSL_SUCCESS;
}

gsl_odeiv_function rhs_func_sin = {
  rhs_sin,
  0
};


int test_rk4(void)
{
  int s = 0;
  double y[2];
  double h = 1.0e-04;
  double t;
  double del;
  double delmax = 0.0;
  int count;

  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_rk4.create(2);

  y[0] = 1.0;
  y[1] = 0.0;

  count = 0;
  for(t=0.0; t<4.0; t += h) {
    gsl_odeiv_step_impl(stepper, t, h, y, &rhs_func_lin);
    if(count % 100 == 0) {
      del = fabs((y[1] - (t+h))/y[1]);
      delmax = GSL_MAX_DBL(del, delmax);
      if( del > (count+1.0) * 2.0e-16 ) {
        printf("LINEAR(%22.18g)  %22.18g  %22.18g  %10.6g\n", t+h, y[1], t+h, del);
	s++;
      }
    }
    count++;
  }

  gsl_odeiv_step_reset(stepper);

  y[0] = 1.0;
  y[1] = 0.0;
  delmax = 0.0;

  count = 0;
  for(t=0.0; t<M_PI; t += h) {
    int stat;
    gsl_odeiv_step_impl(stepper, t, h, y, &rhs_func_sin);
    del = fabs((y[1] - sin(t+h))/y[1]);
    delmax = GSL_MAX_DBL(del, delmax);
    if(count % 100 == 0) {
      if(t < 0.5*M_PI) {
        stat = ( del > (count + 1.0) * 2.0e-16 );
      }
      else if(t < 0.7 * M_PI) {
        stat = ( del > 1.0e-12 );
      }
      else if(t < 0.9 * M_PI) {
        stat = t < ( del > 1.0e-10 );
      }
      else {
        stat = t < ( del > 1.0e-08 );
      }
      if(stat != 0) {
        printf("SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t+h, y[1], sin(t+h), del);
      }
      s += stat;
    }
    count++;
  }
  if(delmax > 1.0e-06) {
    printf("SIN(0 .. M_PI)  delmax = %g\n", delmax);
  }


  for(; t< 100.5 * M_PI; t += h) {
    gsl_odeiv_step_impl(stepper, t, h, y, &rhs_func_sin);
  }
  del = fabs((y[1] - sin(t+h))/y[1]);
  delmax = GSL_MAX_DBL(del, delmax);
  if(del > 1.0e-8) {
    s++;
    printf("SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t+h, y[1], sin(t+h), del);
  }
  if(delmax > 1.0e-06) {
    printf("SIN(0 .. 100.5 M_PI)  delmax = %g\n", delmax);
  }

  return s;
}


int main()
{
  gsl_ieee_env_setup ();

  gsl_test(test_rk4(), "Runge-Kutta 4");

  return gsl_test_summary();
}
