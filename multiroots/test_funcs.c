#include <math.h>
#include <gsl_vector.h>
#include <gsl_matrix.h>

#include "test_funcs.h"

void 
rosenbrock_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, -1.2);
  gsl_vector_set (x, 1, 1.0);
}

int 
rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);

  double y0 = 1 - x0;
  double y1 = 10 * (x1 - x0 * x0);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  
  params = 0 ; /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int 
rosenbrock_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  double x0 = gsl_vector_get (x, 0);

  double df00 = -1;
  double df01 = 0;
  double df10 = -20 * x0;
  double df11 = -10;

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);

  params = 0 ; /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int 
rosenbrock_fdf (const gsl_vector * x, void *params,
		gsl_vector * f, gsl_matrix * df)
{
  rosenbrock_f (x, params, f);
  rosenbrock_df (x, params, df);

  return GSL_SUCCESS;
}
