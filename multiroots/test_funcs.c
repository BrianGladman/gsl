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
  double df11 = 10;

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


void
powellsingular_initpt (gsl_vector * x)
{
  gsl_vector_set (x, 0, 3.0);
  gsl_vector_set (x, 1, -1.0);
  gsl_vector_set (x, 2, 0.0);
  gsl_vector_set (x, 3, 1.0);
}

int 
powellsingular_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double x2 = gsl_vector_get (x, 2);
  double x3 = gsl_vector_get (x, 3);

  double y0 = x0 + 10 * x1;
  double y1 = sqrt(5.0) * (x2 - x3);
  double y2 = pow(x1 - 2 * x2, 2.0);
  double y3 = sqrt(10.0) * pow(x0 - x3, 2.0);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);
  gsl_vector_set (f, 3, y3);
  
  params = 0 ; /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int 
powellsingular_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double x2 = gsl_vector_get (x, 2);
  double x3 = gsl_vector_get (x, 3);

  double df00 = 1, df01 = 10, df02 = 0, df03 = 0;
  double df10 = 0, df11 = 0, df12 = sqrt(5.0), df13 = -df12;
  double df20 = 0, df21 = 2*(x1-2*x2), df22 = -2*df21, df23 = 0;
  double df30 = 2*sqrt(10.0)*(x0-x3), df31 = 0, df32=0, df33 = -df30;

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 0, 2, df02);
  gsl_matrix_set (df, 0, 3, df03);

  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);
  gsl_matrix_set (df, 1, 2, df12);
  gsl_matrix_set (df, 1, 3, df13);

  gsl_matrix_set (df, 2, 0, df20);
  gsl_matrix_set (df, 2, 1, df21);
  gsl_matrix_set (df, 2, 2, df22);
  gsl_matrix_set (df, 2, 3, df23);

  gsl_matrix_set (df, 3, 0, df30);
  gsl_matrix_set (df, 3, 1, df31);
  gsl_matrix_set (df, 3, 2, df32);
  gsl_matrix_set (df, 3, 3, df33);

  params = 0 ; /* avoid warning about unused parameters */

  return GSL_SUCCESS;
}

int 
powellsingular_fdf (const gsl_vector * x, void *params,
		gsl_vector * f, gsl_matrix * df)
{
  powellsingular_f (x, params, f);
  powellsingular_df (x, params, df);

  return GSL_SUCCESS;
}
