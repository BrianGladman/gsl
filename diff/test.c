/* differentiation/test.c
 * 
 * Copyright (C) 2000 David Morrison
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_test.h>

double
f1(double x, void *params)
{
  return exp(x);
}

double
df1(double x, void *params)
{
  return exp(x);
}

double
f2(double x, void *params)
{
  if (x >= 0.0) {
    return x * sqrt(x);
  } else {
    return 0.0;
  }
}

double
df2(double x, void *params)
{
  if (x >= 0.0) {
    return 1.5 * sqrt(x);
  } else {
    return 0.0;
  }
}

double
f3(double x, void *params)
{
  if (x != 0.0) {
    return sin(1/x);
  } else {
    return 0.0;
  }
}

double
df3(double x, void *params)
{
  if (x != 0.0) {
    return -cos(1/x)/(x*x);
  } else {
    return 0.0;
  }
}

double
f4(double x, void *params)
{
  return exp(-x*x);
}

double
df4(double x, void *params)
{
  return -2.0 * x * exp(-x*x);
}

double
f5(double x, void *params)
{
  return x*x;
}

double
df5(double x, void *params)
{
  return 2.0 * x;
}

double
f6(double x, void *params)
{
  return 1.0/x;
}

double
df6(double x, void *params)
{
  return -1.0 / (x * x);
}

double
vf1(const gsl_vector *x, void *params)
{
  return gsl_vector_get(x, 0) +
    pow(gsl_vector_get(x, 1), 2.0) +
    pow(gsl_vector_get(x, 2), 3.0);
}   

void
gvf1(const gsl_vector *x, gsl_vector *g)
{
   gsl_vector_set(g, 0, 1.0);
   gsl_vector_set(g, 1, 2.0 * gsl_vector_get(x,1));
   gsl_vector_set(g, 2, 3.0 * gsl_vector_get(x,2) *
		  gsl_vector_get(x,2)); 
}   

double
vf2(const gsl_vector *x, void *params)
{
  return 1.0 * exp(gsl_vector_get(x,0)) +
    2.0 * exp(gsl_vector_get(x,1)) +
    3.0 * exp(gsl_vector_get(x,2));
}   

void
gvf2(const gsl_vector *x, gsl_vector *g)
{
   gsl_vector_set(g, 0, 1.0 * exp(gsl_vector_get(x,0)));
   gsl_vector_set(g, 1, 2.0 * exp(gsl_vector_get(x,1)));
   gsl_vector_set(g, 2, 3.0 * exp(gsl_vector_get(x,2)));
}   

int
main ()
{ 
  int i;
  double x, exp_result;
  double result, abserr;
  gsl_vector *vx, *vg;
  gsl_multimin_function vf = {0,0,0};
  gsl_function f = {0,0}, df = {0,0}; 

  f.function = &f1;
  df.function = &df1;
  x = 1.0;
  exp_result = GSL_FN_EVAL(&df, x);
  gsl_diff_central(&f, x, &result, &abserr);
  printf("result = %12.8G\n", result);
  printf("abserr = %12.8G\n", abserr);
  gsl_test_rel(result, exp_result, abserr,
	       "central difference method");

  f.function = &f3;
  df.function = &df3;
  x = 0.45;
  exp_result = GSL_FN_EVAL(&df, x);
  gsl_diff_central(&f, x, &result, &abserr);
  printf("result = %12.8G\n", result);
  printf("abserr = %12.8G\n", abserr);
  gsl_test_rel(result, exp_result, abserr,
	       "central difference method");

  f.function = &f4;
  df.function = &df4;
  x = 0.5;
  exp_result = GSL_FN_EVAL(&df, x);
  gsl_diff_central(&f, x, &result, &abserr);
  printf("result = %12.8G\n", result);
  printf("abserr = %12.8G\n", abserr);
  gsl_test_rel(result, exp_result, abserr,
	       "central difference method");

  f.function = &f2;
  df.function = &df2;
  x = 0.1;
  exp_result = GSL_FN_EVAL(&df, x);
  gsl_diff_forward(&f, x, &result, &abserr);
  printf("result = %12.8G\n", result);
  printf("abserr = %12.8G\n", abserr);
  gsl_test_rel(result, exp_result, abserr,
	       "forward difference method");

  f.function = &f5;
  df.function = &df5;
  x = 0.0;
  exp_result = GSL_FN_EVAL(&df, x);
  gsl_diff_forward(&f, x, &result, &abserr);
  printf("result = %12.8G\n", result);
  printf("abserr = %12.8G\n", abserr);
  gsl_test_rel(result, exp_result, abserr,
	       "forward difference method");

  f.function = &f6;
  df.function = &df6;
  x = 10.0;
  exp_result = GSL_FN_EVAL(&df, x);
  gsl_diff_backward(&f, x, &result, &abserr);
  printf("result = %12.8G\n", result);
  printf("abserr = %12.8G\n", abserr);
  gsl_test_rel(result, exp_result, abserr,
	       "forward backward method");

  vf.f = &vf1;
  vf.n = 3;
  vx = gsl_vector_calloc(3);
  for (i = 0; i < 3; i++) 
    {
      gsl_vector_set(vx, i, 1.0);
    }
  vg = gsl_vector_calloc(3);
  gsl_diff_gradient(&vf, vx, vg);
  printf("calc g = (%12.8G,%12.8G,%12.8G)\n", 
	 gsl_vector_get(vg, 0),
	 gsl_vector_get(vg, 1),
	 gsl_vector_get(vg, 2));
  gvf1(vx, vg);
  printf("true g = (%12.8G,%12.8G,%12.8G)\n", 
	 gsl_vector_get(vg, 0),
	 gsl_vector_get(vg, 1),
	 gsl_vector_get(vg, 2));
  
  vf.f = &vf2;
  vf.n = 3;
  for (i = 0; i < 3; i++) 
    {
      gsl_vector_set(vx, i, i);
    }
  gsl_diff_gradient(&vf, vx, vg);
  printf("calc g = (%12.8G,%12.8G,%12.8G)\n", 
	 gsl_vector_get(vg, 0),
	 gsl_vector_get(vg, 1),
	 gsl_vector_get(vg, 2));
  gvf2(vx, vg);
  printf("true g = (%12.8G,%12.8G,%12.8G)\n", 
	 gsl_vector_get(vg, 0),
	 gsl_vector_get(vg, 1),
	 gsl_vector_get(vg, 2));
  
  return gsl_test_summary();
}
