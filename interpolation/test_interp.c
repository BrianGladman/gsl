/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl_test.h>
#include <gsl_errno.h>
#include <gsl_interp.h>


double test_func_1(double x) { return sin(x); }
double test_func_2(double x) { return x*x; }


typedef  struct _xy_table  xy_table;

struct _xy_table {
  double * x;
  double * y;
  int      n;
};


int
test_interp(xy_table * tab, gsl_interp_factory * fact, double (*f)(double))
{
  gsl_interp_accel * a = gsl_interp_accel_new();
  gsl_interp_obj * interp = fact->create(tab->x, tab->y, tab->n);
  
  
  gsl_interp_accel_free(a);
  gsl_interp_obj_free(interp);
  
  return 0;
}


int main(int argc, char ** argv)
{
  int status = 0;
  int s;
  int i;

  double x;
  double xmin = 0.0;
  double xmax = 2.0*M_PI;
  double dx = 1.e-4 * xmax;

  double (*test_func)(double) = test_func_1;

  xy_table data_table;
  data_table.n = 200;
  data_table.x = (double *) malloc(data_table.n * sizeof(double));
  data_table.y = (double *) malloc(data_table.n * sizeof(double));

  for(i=0; i<data_table.n; i++) {
    data_table.x[i] = xmin + (xmax - xmin) * i/(double)(data_table.n-1);
    data_table.y[i] = test_func(data_table.x[i]);
  }

  s = 0;
  s += test_interp(&data_table, &gsl_interp_factory_linear, test_func);
  gsl_test(s, "linear interpolation");
  status += s;

  s = 0;
  s += test_interp(&data_table, &gsl_interp_factory_cspline_natural, test_func);
  gsl_test(s, "cspline interpolation");
  status += s;

  return gsl_test_summary();
}
