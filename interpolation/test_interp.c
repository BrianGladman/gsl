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


typedef  double TEST_FUNC(double);
typedef  struct _xy_table    xy_table;


struct _xy_table {
  double * x;
  double * y;
  size_t n;
};


static void
alloc_xy_table(xy_table * data_table, size_t size)
{
  data_table->n = size;
  data_table->x = (double *) malloc(data_table->n * sizeof(double));
  data_table->y = (double *) malloc(data_table->n * sizeof(double));
}


static int
test_interp(xy_table * data_table, const gsl_interp_factory * factory, xy_table * test_table)
{
  int status = 0;
  size_t i;

  gsl_interp_accel * a      = gsl_interp_accel_new();
  gsl_interp_obj   * interp = factory->create(data_table->x,
                                              data_table->y,
					      data_table->n);
  
  for(i=0; i<test_table->n; i++) {
    double x = test_table->x[i];
    double y;
    double diff;
    gsl_interp_eval_impl(interp, data_table->x, data_table->y, x, a, &y);
    diff = y - test_table->y[i];
    if(fabs(diff) > 1.e-10) status++;
  }
  
  gsl_interp_accel_free(a);
  gsl_interp_obj_free(interp);
  
  return status;
}


static int 
test_linear(void)
{
  int s;
  
  double data_x[3] = { 0.0, 1.0, 2.0 };
  double data_y[3] = { 0.0, 1.0, 2.0 };
  double test_x[1] = { 0.0 };
  double test_y[1] = { 0.0 };

  xy_table  data_table = { data_x, data_y, sizeof(data_x)/sizeof(double) };
  xy_table  test_table = { test_x, test_y, sizeof(test_x)/sizeof(double) };
  
  s = test_interp(&data_table, &gsl_interp_factory_linear, &test_table);
  gsl_test(s, "linear interpolation");
  return s;
}

static int 
test_cspline_natural(void)
{ 
  int s ;
  
  double data_x[3] = { 0.0, 1.0, 2.0 };
  double data_y[3] = { 0.0, 1.0, 2.0 };
  double test_x[1] = { 0.0 };
  double test_y[1] = { 0.0 };

  xy_table  data_table = { data_x, data_y, sizeof(data_x)/sizeof(double) };
  xy_table  test_table = { test_x, test_y, sizeof(test_x)/sizeof(double) };
  
  s = test_interp(&data_table, &gsl_interp_factory_cspline_natural, &test_table);
  gsl_test(s, "cspline interpolation");
  return s;
}


int main(int argc, char ** argv)
{
  int status = 0;
  
  argc = 0;  /* prevent warnings about unused parameters */
  argv = 0;

  status += test_linear();
  status += test_cspline_natural();

  return gsl_test_summary();
}



