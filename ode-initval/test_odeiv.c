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

gsl_odeiv_system rhs_func_lin = {
  rhs_linear,
  0,
  0
};


/* RHS for sin(t),cos(t) */

int rhs_sin(double t, const double y[], double f[], void * params)
{
  f[0] = -y[1];
  f[1] =  y[0];
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_sin = {
  rhs_sin,
  0,
  0
};


/* RHS for a exp(t)+ b exp(-t) */

int rhs_exp(double t, const double y[], double f[], void * params)
{
  f[0] = y[1];
  f[1] = y[0];
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_exp = {
  rhs_exp,
  0,
  0
};


/* RHS for stiff example */

int rhs_stiff(double t, const double y[], double f[], void * params)
{
  f[0] =  998.0*y[0] + 1998.0*y[1];
  f[1] = -999.0*y[0] - 1999.0*y[1];
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_stiff = {
  rhs_stiff,
  0,
  0
};


int test_stepper_linear(gsl_odeiv_step * stepper, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del;
  double delmax = 0.0;
  int count = 0;

  y[0] = 1.0;
  y[1] = 0.0;

  for(t=0.0; t<4.0; t += h) {
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, &rhs_func_lin);
    if(count % 100 == 0) {
      del = fabs((y[1] - (t+h))/y[1]);
      delmax = GSL_MAX_DBL(del, delmax);
      if(del > (count+1.0) * base_prec) {
        printf("  LINEAR(%20.17g)  %20.17g  %20.17g  %8.4g\n", t+h, y[1], t+h, del);
	s++;
      }
    }
    count++;
  }

  return s;
}


int test_stepper_sin(gsl_odeiv_step * stepper, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del;
  double delmax = 0.0;
  int count = 0;

  y[0] = 1.0;
  y[1] = 0.0;

  for(t=0.0; t<M_PI; t += h) {
    int stat;
    double sin_th = sin(t+h);
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, &rhs_func_sin);
    del = fabs((y[1] - sin_th)/y[1]);
    delmax = GSL_MAX_DBL(del, delmax);
    {
      if(t < 0.5*M_PI) {
        stat = ( del > (count + 1.0) * base_prec );
      }
      else if(t < 0.7 * M_PI) {
        stat = ( del > 1.0e+04 * base_prec );
      }
      else if(t < 0.9 * M_PI) {
        stat = ( del > 1.0e+06 * base_prec );
      }
      else {
        stat = ( del > 1.0e+09 * base_prec );
      }
      if(stat != 0) {
        printf("  SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t+h, y[1], sin_th, del);
      }
      s += stat;
    }
    count++;
  }
  if(delmax > 1.0e+09 * base_prec) {
    s++;
    printf("  SIN(0 .. M_PI)  delmax = %g\n", delmax);
  }


  for(; t< 100.5 * M_PI; t += h) {
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, &rhs_func_sin);
    count++;
  }
  del = fabs((y[1] - sin(t))/y[1]);
  if(del > count * 2.0 * base_prec) {
    s++;
    printf("  SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t+h, y[1], sin(t), del);
  }

  return s;
}


int test_stepper_exp(gsl_odeiv_step * stepper, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del;
  int count = 0;

  y[0] = 1.0;
  y[1] = 1.0;

  for(t=0.0; t<20.0; t += h) {
    double ex = exp(t + h);
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, &rhs_func_exp);
    del = fabs((y[1] - ex)/y[1]);
    if(del > (count+1.0) * 2.0 * base_prec) {
      printf("  EXP(%20.17g)  %20.17g  %20.17g  %8.4g\n", t+h, y[1], ex, del);
      s++;
    }
    count++;
  }

  return s;
}


int test_stepper_stiff(gsl_odeiv_step * stepper, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del;
  int count = 0;

  y[0] = 1.0;
  y[1] = 0.0;

  for(t=0.0; t<20.0; t += h) {
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, &rhs_func_stiff);
    if(t > 0.04) {
      double arg = t + h;
      double e1 = exp(-arg);
      double e2 = exp(-1000.0*arg);
      double u = 2.0*e1 - e2;
      /* double v = -e1 + e2; */
      del = fabs((y[0] - u)/y[0]);
      if(del > (count+1.0) * 100.0 * base_prec) {
        printf("  STIFF(%20.17g)  %20.17g  %20.17g  %8.4g\n", arg, y[0], u, del);
	s++;
      }
    }
    count++;
  }

  return s;
}


int test_stepper_rk2(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_rk2.create(2);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-03, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 4.0e-04, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;

  s = test_stepper_exp(stepper, 4.0e-04, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;

  gsl_odeiv_step_free(stepper);
  return stat;
}


int test_stepper_rk4(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_rk4.create(2);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-03, GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 2.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;

  s = test_stepper_exp(stepper, 2.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;

  gsl_odeiv_step_free(stepper);
  return stat;
}


int test_stepper_rkck(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_rkck.create(2);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-03, GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 5.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 5.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;

  gsl_odeiv_step_free(stepper);
  return stat;
}


int test_stepper_rk8pd(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_rk8pd.create(2);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 5.0e-02, GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 5.0e-02, GSL_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 5.0e-02, GSL_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);


  s = test_stepper_sin(stepper, 0.2, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 0.2, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);


  gsl_odeiv_step_free(stepper);
  return stat;
}


int test_stepper_rk4imp(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_rk4imp.create(2);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-03, GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 2.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 2.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 2.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;

  gsl_odeiv_step_free(stepper);
  return stat;
}


int test_stepper_gear1(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_gear1.create(2);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-03, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 4.0e-04, 10.0 * GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 1.0e-03, 20.0 * GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 1.0e-03, 10.0 * GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;

  gsl_odeiv_step_free(stepper);
  return stat;
}


int test_stepper_gear2(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_factory_gear2.create(2);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 2.0e-04, 20.0 * GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 2.0e-04, 20.0 * GSL_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 2.0e-04, 20.0 * GSL_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 2.0e-04, 20.0 * GSL_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;

  gsl_odeiv_step_free(stepper);
  return stat;
}


int main()
{
  gsl_ieee_env_setup();

  gsl_test(test_stepper_rk2(),     "Runge-Kutta 2(3), Euler-Cauchy");
  gsl_test(test_stepper_rk4(),     "Runge-Kutta 4, Classical");
  gsl_test(test_stepper_rkck(),    "Runge-Kutta 4(5), Cash-Karp");
  gsl_test(test_stepper_rk8pd(),   "Runge-Kutta 8(9), Prince-Dormand");
  gsl_test(test_stepper_rk4imp(),  "Runge-Kutta 4, Gaussian implicit");
  gsl_test(test_stepper_gear1(),   "Gear 1");
/*  gsl_test(test_stepper_gear2(),   "Gear 2"); */

  return gsl_test_summary();
}
