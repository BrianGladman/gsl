#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

int
brown_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double x2 = gsl_vector_get (x, 2);
  double x3 = gsl_vector_get (x, 3);
  size_t i;

  for (i = 0; i < 20; i++)
    {
      double ti = 0.2 * (i+1);
      double ui = x0 + x1 * ti - exp(ti);
      double vi = x2 + x3 * sin(ti) - cos(ti);

      gsl_vector_set (f, i, ui*ui + vi*vi);      
    }

  return GSL_SUCCESS;
}

int
brown_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);
  double x2 = gsl_vector_get (x, 2);
  double x3 = gsl_vector_get (x, 3);
  size_t i;

  for (i = 0; i < 20; i++)
    {
      double ti = 0.2 * (i+1);
      double ui = x0 + x1 * ti - exp(ti);
      double vi = x2 + x3 * sin(ti) - cos(ti);
      
      gsl_matrix_set (df, i, 0, 2*ui);
      gsl_matrix_set (df, i, 1, 2*ui * ti);
      gsl_matrix_set (df, i, 2, 2*vi);
      gsl_matrix_set (df, i, 3, 2*vi * sin(ti));

    }
  return GSL_SUCCESS;
}

int
brown_fdf (const gsl_vector * x, void *params,
		gsl_vector * f, gsl_matrix * df)
{
  brown_f (x, params, f);
  brown_df (x, params, df);

  return GSL_SUCCESS;
}

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 20;
  const size_t p = 4;

  gsl_multifit_function_fdf f = {&brown_f, &brown_df, &brown_fdf, n, p, 0};

  double x_init[4] = {25, 5, -5, -1};
  gsl_vector x = gsl_vector_view (x_init, p);

  T = gsl_multifit_fdfsolver_lmder;
  s = gsl_multifit_fdfsolver_alloc (T, &f, &x);

//  gsl_multifit_fdfsolver_set (

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      print_state (iter, s);

      //if (status)
	//break;

      //status = gsl_multifit_test_residual (s->f, 0.0000001);
    }
  while (iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);

}

int
print_state (size_t iter,   gsl_multifit_fdfsolver * s)
{
  printf ("iter = %3u x = % 15.8f % 15.8f % 15.8f % 15.8f |f(x)| = %g\n",
	  iter,
	  gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2), gsl_vector_get (s->x, 3),
          gsl_blas_dnrm2 (s->f));
}


