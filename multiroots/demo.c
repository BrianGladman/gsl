#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct rparams
  {
    double a;
    double b;
  };

int
rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double a = ((struct rparams *) params)->a;
  double b = ((struct rparams *) params)->b;

  double x0 = gsl_vector_get (x, 0);
  double x1 = gsl_vector_get (x, 1);

  double y0 = a * (1 - x0);
  double y1 = b * (x1 - x0 * x0);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}

#ifdef DERIV
int
rosenbrock_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  double a = ((struct rparams *) params)->a;
  double b = ((struct rparams *) params)->b;

  double x0 = gsl_vector_get (x, 0);

  double df00 = -a;
  double df01 = 0;
  double df10 = -2 * b  * x0;
  double df11 = b;

  gsl_matrix_set (df, 0, 0, df00);
  gsl_matrix_set (df, 0, 1, df01);
  gsl_matrix_set (df, 1, 0, df10);
  gsl_matrix_set (df, 1, 1, df11);

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
#endif

#ifdef DERIV
#define SOLVER gsl_multiroot_fdfsolver
#define SOLVER_TYPE gsl_multiroot_fdfsolver_type
#else
#define SOLVER gsl_multiroot_fsolver
#define SOLVER_TYPE gsl_multiroot_fsolver_type
#endif

int
main (void)
{
  const SOLVER_TYPE *T;
  SOLVER *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 2;
  struct rparams p = {1.0, 10.0};
#ifdef DERIV
  gsl_multiroot_function_fdf f = {&rosenbrock_f, &rosenbrock_df, &rosenbrock_fdf, n, &p};
#else
  gsl_multiroot_function f = {&rosenbrock_f, n, &p};
#endif

  double x_init[2] = {-10.0, -5.0};
  gsl_vector *x = gsl_vector_alloc (n);

  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

#ifdef DERIV
  T = gsl_multiroot_fdfsolver_gnewton;
  s = gsl_multiroot_fdfsolver_alloc (T, &f, x);
#else
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, &f, x);
#endif

  print_state (iter, s);

  do
    {
      iter++;
#ifdef DERIV
      status = gsl_multiroot_fdfsolver_iterate (s);
#else
      status = gsl_multiroot_fsolver_iterate (s);
#endif

      print_state (iter, s);

      if (status)
        break;

      status = gsl_multiroot_test_residual (s->f, 0.0000001);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

#ifdef DERIV
  gsl_multiroot_fdfsolver_free (s);
#else
  gsl_multiroot_fsolver_free (s);
#endif

  gsl_vector_free (x);
}

int
print_state (size_t iter, SOLVER * s)
{
  printf ("iter = %3u x = % 15.8f % 15.8f  f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
}


