/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl_test.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_matrix.h>
#include <gsl_linalg.h>
#include <gsl_ieee_utils.h>
#include "gsl_odeiv.h"


/* RHS for a + b t */

int rhs_linear(double t, const double y[], double f[], void * params)
{
  f[0] = 0.0;
  f[1] = y[0];
  return GSL_SUCCESS;
}

int jac_linear(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.dim2  = 2;
  dfdy_mat.data  = dfdy;
  dfdy_mat.block = 0;
  gsl_matrix_set(&dfdy_mat, 0, 0, 0.0);
  gsl_matrix_set(&dfdy_mat, 0, 1, 0.0);
  gsl_matrix_set(&dfdy_mat, 1, 0, 1.0);
  gsl_matrix_set(&dfdy_mat, 1, 1, 0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_lin = {
  rhs_linear,
  jac_linear,
  2,
  0
};


/* RHS for sin(t),cos(t) */

int rhs_sin(double t, const double y[], double f[], void * params)
{
  f[0] = -y[1];
  f[1] =  y[0];
  return GSL_SUCCESS;
}

int jac_sin(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.data  = dfdy;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.dim2  = 2;
  dfdy_mat.block = 0;
  gsl_matrix_set(&dfdy_mat, 0, 0,  0.0);
  gsl_matrix_set(&dfdy_mat, 0, 1, -1.0);
  gsl_matrix_set(&dfdy_mat, 1, 0,  1.0);
  gsl_matrix_set(&dfdy_mat, 1, 1,  0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_sin = {
  rhs_sin,
  0,
  2,
  0
};


/* RHS for a exp(t)+ b exp(-t) */

int rhs_exp(double t, const double y[], double f[], void * params)
{
  f[0] = y[1];
  f[1] = y[0];
  return GSL_SUCCESS;
}

int jac_exp(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.data = dfdy;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.dim2  = 2;
  dfdy_mat.block = 0;
  gsl_matrix_set(&dfdy_mat, 0, 0, 0.0);
  gsl_matrix_set(&dfdy_mat, 0, 1, 1.0);
  gsl_matrix_set(&dfdy_mat, 1, 0, 1.0);
  gsl_matrix_set(&dfdy_mat, 1, 1, 0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_exp = {
  rhs_exp,
  jac_exp,
  2,
  0
};


/* RHS for stiff example */

int rhs_stiff(double t, const double y[], double f[], void * params)
{
  f[0] =  998.0*y[0] + 1998.0*y[1];
  f[1] = -999.0*y[0] - 1999.0*y[1];
  return GSL_SUCCESS;
}

int jac_stiff(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.data = dfdy;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.dim2  = 2;
  dfdy_mat.block = 0;
  gsl_matrix_set(&dfdy_mat, 0, 0,   998.0);
  gsl_matrix_set(&dfdy_mat, 0, 1,  1998.0);
  gsl_matrix_set(&dfdy_mat, 1, 0,  -999.0);
  gsl_matrix_set(&dfdy_mat, 1, 1, -1999.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_stiff = {
  rhs_stiff,
  jac_stiff,
  2,
  0
};


int test_stepper_linear(gsl_odeiv_step * stepper, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del;
  int count = 0;

  y[0] = 1.0;
  y[1] = 0.0;

  for(t=0.0; t<4.0; t += h) {
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, 0, 0, &rhs_func_lin);
    /* if(count % 100 == 0) */ {
      del = fabs((y[1] - (t+h))/y[1]);
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
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, 0, 0, &rhs_func_sin);
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
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, 0, 0, &rhs_func_sin);
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
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, 0, 0, &rhs_func_exp);
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
    gsl_odeiv_step_impl(stepper, t, h, y, yerr, 0, 0, &rhs_func_stiff);
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
  gsl_odeiv_step * stepper = gsl_odeiv_step_rk2_new();
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
  gsl_odeiv_step * stepper = gsl_odeiv_step_rk4_new();
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
  gsl_odeiv_step * stepper = gsl_odeiv_step_rkck_new();
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
  gsl_odeiv_step * stepper = gsl_odeiv_step_rk8pd_new();
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


  s = test_stepper_sin(stepper, 0.25, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 0.25, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);


  gsl_odeiv_step_free(stepper);

  return stat;
}


int test_stepper_rk2imp(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_rk2imp_new();
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-02, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 3.0e-04, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 5.0e-03, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 1.0e-03, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;

  gsl_odeiv_step_free(stepper);

  return stat;
}


int test_stepper_rk4imp(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_rk4imp_new();
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-02, GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 4.0e-04, GSL_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 1.0e-03, GSL_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 1.0e-03, GSL_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;

  gsl_odeiv_step_free(stepper);

  return stat;
}


int test_stepper_gear1(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_gear1_new();
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
  gsl_odeiv_step * stepper = gsl_odeiv_step_gear2_new();
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-03, GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_sin(stepper, 4.0e-04, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  SIN");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_exp(stepper, 2.0e-04,  1.0e+05 * GSL_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 2.0e-04, 1.0e+05 * GSL_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;


  s = test_stepper_exp(stepper, 3.0e-03,  GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 1.0e-03, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;


  gsl_odeiv_step_free(stepper);

  return stat;
}


int test_stepper_bsimp(void)
{
  gsl_odeiv_step * stepper = gsl_odeiv_step_bsimp_new(1.0e-04);
  int stat = 0;
  int s;

  s = test_stepper_linear(stepper, 1.0e-03, 1.0e+04 * GSL_DBL_EPSILON);
  gsl_test(s, "  LINEAR");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 2.0e-04, 1.0e+05 * GSL_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;

  s = test_stepper_exp(stepper, 3.0e-03,  GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  EXP");
  stat += s;
  gsl_odeiv_step_reset(stepper);

  s = test_stepper_stiff(stepper, 1.0e-03, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  STIFF");
  stat += s;


  gsl_odeiv_step_free(stepper);

  return stat;
}


int test_evolve_system_flat(
  gsl_odeiv_step * step, 
  gsl_odeiv_evolve_mon * mon,
  const gsl_odeiv_system * sys,
  double t0, double t1, double hstart,
  double y[], double yfin[],
  double err_target)
{
  int s = 0;
  double frac;

  gsl_odeiv_evolve * e = gsl_odeiv_evolve_new();

  gsl_odeiv_evolve_impl(e, mon, 0, step, sys, t0, t1, hstart, y);

  frac = fabs((y[1] - yfin[1])/y[1]) + fabs((y[0] - yfin[0])/y[0]);
  if(frac > 2.0 * e->count * err_target) {
    s++;
  }

  gsl_odeiv_evolve_free(e);

  return s;
}


int test_evolve_system(
  gsl_odeiv_step * step, 
  gsl_odeiv_evolve_mon * mon,
  const gsl_odeiv_system * sys,
  double t0, double t1, double hstart,
  double y[], double yfin[],
  double err_target)
{
  int s = 0;
  double frac;

  gsl_odeiv_evolve_control * c = gsl_odeiv_evolve_control_standard_new(err_target, 0.0, 1.0, 1.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_new();

  gsl_odeiv_evolve_impl(e, mon, c, step, sys, t0, t1, hstart, y);

  frac = fabs((y[1] - yfin[1])/y[1]) + fabs((y[0] - yfin[0])/y[0]);
  if(frac > 2.0 * e->count * err_target) {
    s++;
  }

/* printf(" COUNT= %d    STUT= %d\n", e->count, e->count_stutter); */

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_evolve_control_free(c);
  
  return s;
}


int test_evolve(void)
{
  int s;
  int status = 0;
  double y[2];
  double yfin[2];

  gsl_odeiv_evolve_mon * mon = 0 /* gsl_odeiv_evolve_mon_stream_new(stdout) */;
  

  gsl_odeiv_step * step = gsl_odeiv_step_rkck_new();
  y[0] = 1.0;
  y[1] = 1.0;
  yfin[0] = exp(10.0);
  yfin[1] = yfin[0];
  s = test_evolve_system(step, mon, &rhs_func_exp, 0.0, 10.0, 1.0e-03, y, yfin, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  Evolve: exp, rkck");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_rk8pd_new();
  y[0] = 1.0;
  y[1] = 1.0;
  yfin[0] = exp(10.0);
  yfin[1] = yfin[0];
  s = test_evolve_system(step, mon, &rhs_func_exp, 0.0, 10.0, 1.0e-03, y, yfin, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  Evolve: exp, rk8pd");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_rk8pd_new();
  y[0] = 1.0;
  y[1] = 0.0;
  yfin[0] = cos(2.0);
  yfin[1] = sin(2.0);
  s = test_evolve_system(step, mon, &rhs_func_sin, 0.0, 2.0, 1.0e-03, y, yfin, GSL_SQRT_DBL_EPSILON);
  gsl_test(s, "  Evolve: sin, rk8pd");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_rk2imp_new();
  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 5.0;
    double e1 = exp(-arg);
    double e2 = exp(-1000.0*arg);
    yfin[0] = 2.0*e1 - e2;
    yfin[1] = -e1 + e2;
  }
  s = test_evolve_system(step, mon, &rhs_func_stiff, 0.0, 5.0, 2.0e-03, y, yfin, 1.0e-05);
  gsl_test(s, "  Evolve: stiff, rk2imp");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_rk4imp_new();
  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 5.0;
    double e1 = exp(-arg);
    double e2 = exp(-1000.0*arg);
    yfin[0] = 2.0*e1 - e2;
    yfin[1] = -e1 + e2;
  }
  s = test_evolve_system(step, mon, &rhs_func_stiff, 0.0, 5.0, 2.0e-03, y, yfin, 1.0e-05);
  gsl_test(s, "  Evolve: stiff, rk4imp");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_gear1_new();
  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 1.0;
    double e1 = exp(-arg);
    double e2 = exp(-1000.0*arg);
    yfin[0] = 2.0*e1 - e2;
    yfin[1] = -e1 + e2;
  }
  s = test_evolve_system_flat(step, mon, &rhs_func_stiff, 0.0, 1.0, 5.0e-03, y, yfin, 1.0e-04);
  gsl_test(s, "  Evolve Flat: stiff, gear1");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_gear1_new();
  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 5.0;
    double e1 = exp(-arg);
    double e2 = exp(-1000.0*arg);
    yfin[0] = 2.0*e1 - e2;
    yfin[1] = -e1 + e2;
  }
  s = test_evolve_system(step, mon, &rhs_func_stiff, 0.0, 5.0, 5.0e-03, y, yfin, 1.0e-04);
  gsl_test(s, "  Evolve: stiff, gear1");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_gear2_new();
  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 1.0;
    double e1 = exp(-arg);
    double e2 = exp(-1000.0*arg);
    yfin[0] = 2.0*e1 - e2;
    yfin[1] = -e1 + e2;
  }
  s = test_evolve_system_flat(step, mon, &rhs_func_stiff, 0.0, 1.0, 5.0e-04, y, yfin, 1.0e-06);
  gsl_test(s, "  Evolve Flat: stiff, gear2");
  status += s;
  gsl_odeiv_step_free(step);


  step = gsl_odeiv_step_gear2_new();
  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 5.0;
    double e1 = exp(-arg);
    double e2 = exp(-1000.0*arg);
    yfin[0] = 2.0*e1 - e2;
    yfin[1] = -e1 + e2;
  }
  s = test_evolve_system(step, mon, &rhs_func_stiff, 0.0, 5.0, 2.0e-03, y, yfin, 1.0e-05);
  gsl_test(s, "  Evolve: stiff, gear2");
  status += s;
  gsl_odeiv_step_free(step);


  gsl_odeiv_evolve_mon_free(mon);

  return status;
}



int main()
{
  gsl_ieee_env_setup();

  gsl_test(test_stepper_rk2(),     "Runge-Kutta 2(3), Euler-Cauchy");
  gsl_test(test_stepper_rk4(),     "Runge-Kutta 4, Classical");
  gsl_test(test_stepper_rkck(),    "Runge-Kutta 4(5), Cash-Karp");
  gsl_test(test_stepper_rk8pd(),   "Runge-Kutta 8(9), Prince-Dormand");
  gsl_test(test_stepper_rk2imp(),  "Runge-Kutta 2, Gaussian implicit");
  gsl_test(test_stepper_rk4imp(),  "Runge-Kutta 4, Gaussian implicit");
  gsl_test(test_stepper_gear1(),   "Gear 1");
  gsl_test(test_stepper_gear2(),   "Gear 2");

  gsl_test(test_stepper_bsimp(),   "Bulirsch-Stoer Implicit");


  gsl_test(test_evolve(),  "Evolution");

  return gsl_test_summary();
}
