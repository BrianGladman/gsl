
/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

typedef double TEST_FUNC (double);
typedef struct _xy_table xy_table;

struct _xy_table
  {
    double *x;
    double *y;
    size_t n;
  };

static int
test_interp (
  xy_table * data_table,
  const gsl_interp_factory * factory,
  xy_table * test_table,
  xy_table * test_d_table,
  xy_table * test_i_table)
{
  int status = 0;
  size_t i;

  gsl_interp_accel *a = gsl_interp_accel_new ();
  gsl_interp_obj *interp = factory->create (data_table->x,
					    data_table->y,
					    data_table->n);

  for (i = 0; i < test_table->n; i++)
    {
      double x = test_table->x[i];
      double y;
      double deriv;
      double integ;
      double diff_y, diff_deriv, diff_integ;
      gsl_interp_eval_impl (interp, data_table->x, data_table->y, x, a, &y);
      gsl_interp_eval_deriv_impl (interp, data_table->x, data_table->y, x, a, &deriv);
      gsl_interp_eval_integ_impl (interp, data_table->x, data_table->y, 0.0, x, a, &integ);
      diff_y = y - test_table->y[i];
      diff_deriv = deriv - test_d_table->y[i];
      diff_integ = integ - test_i_table->y[i];
      if (fabs (diff_y) > 1.e-10 || fabs(diff_deriv) > 1.0e-10 || fabs(diff_integ) > 1.0e-10) {
	status++;
      }
    }

  gsl_interp_accel_free (a);
  gsl_interp_obj_free (interp);

  return status;
}


static int
test_linear (void)
{
  int s;

  double data_x[4] =
  {0.0, 1.0, 2.0, 3.0};
  double data_y[4] =
  {0.0, 1.0, 2.0, 3.0};
  double test_x[4] =
  {0.0, 0.5, 1.5, 2.5};
  double test_y[4] =
  {0.0, 0.5, 1.5, 2.5};
  double test_dy[4] =
  {1.0, 1.0, 1.0, 1.0};
  double test_iy[4] =
  {0.0, 0.125, 9.0/8.0, 25.0/8.0};

  xy_table data_table =
  {data_x, data_y, sizeof (data_x) / sizeof (double)};
  xy_table test_table =
  {test_x, test_y, sizeof (test_x) / sizeof (double)};
  xy_table test_d_table =
  {test_x, test_dy, sizeof (test_x) / sizeof (double)};
  xy_table test_i_table =
  {test_x, test_iy, sizeof (test_x) / sizeof (double)};

  s = test_interp (&data_table, &gsl_interp_factory_linear, &test_table, &test_d_table, &test_i_table);
  gsl_test (s, "linear interpolation");
  return s;
}

static int
test_cspline_natural (void)
{
  int s;

  double data_x[3] =
  {0.0, 1.0, 2.0};
  double data_y[3] =
  {0.0, 1.0, 2.0};
  double test_x[2] =
  {0.0, 0.5};
  double test_y[2] =
  {0.0, 0.5};
  double test_dy[2] =
  {1.0, 1.0};
  double test_iy[2] =
  {0.0, 0.125};

  xy_table data_table =
  {data_x, data_y, sizeof (data_x) / sizeof (double)};
  xy_table test_table =
  {test_x, test_y, sizeof (test_x) / sizeof (double)};
  xy_table test_d_table =
  {test_x, test_dy, sizeof (test_x) / sizeof (double)};
  xy_table test_i_table =
  {test_x, test_iy, sizeof (test_x) / sizeof (double)};
  
  s = test_interp (&data_table, &gsl_interp_factory_cspline_natural, &test_table, &test_d_table, &test_i_table);
  gsl_test (s, "cspline interpolation");
  return s;
}


int 
main (int argc, char **argv)
{
  int status = 0;

  argc = 0;			/* prevent warnings about unused parameters */
  argv = 0;

  status += test_linear ();
  status += test_cspline_natural ();

  return gsl_test_summary ();
}
