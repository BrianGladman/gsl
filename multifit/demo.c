#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

struct data {
  size_t n;
  double * y;
  double * sigma;
};

int
expb_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  size_t n = ((struct data *)params)->n;
  double *y = ((struct data *)params)->y;
  double *sigma = ((struct data *) params)->sigma;

  double x0 = gsl_vector_get (x, 0);
  double A = gsl_vector_get (x, 1);
  double lambda = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      double t = i;
      double Yi = x0 + A * exp (lambda * t);
      gsl_vector_set (f, i, (y[i] - Yi)/sigma[i]);
    }

  return GSL_SUCCESS;
}

int
expb_df (const gsl_vector * x, void *params, gsl_matrix * df)
{
  size_t n = ((struct data *)params)->n;
  double *sigma = ((struct data *) params)->sigma;

  double A = gsl_vector_get (x, 1);
  double lambda = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++)
    {
      /*  Yi = x0 + A * exp(lambda * i)  */
      double t = i;
      double s = sigma[i];
      gsl_matrix_set (df, i, 0, -1/s);
      gsl_matrix_set (df, i, 1, -exp (lambda * t)/s);
      gsl_matrix_set (df, i, 2, -t * A * exp (lambda * t)/s);

    }
  return GSL_SUCCESS;
}

int
expb_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * df)
{
  expb_f (x, params, f);
  expb_df (x, params, df);

  return GSL_SUCCESS;
}

int
main (void)
{
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 40;
  const size_t p = 3;

  gsl_vector *g = gsl_vector_calloc (p);
  gsl_matrix *covar = gsl_matrix_alloc (p, p);

  double y[n], sigma[n];

  struct data d = { n, y, sigma};
  
  gsl_multifit_function_fdf f =
    { &expb_f, &expb_df, &expb_fdf, n, p, (void *) &d };

  double x_init[3] = { 0.0, 1.0, -0.0 };

  gsl_vector x = gsl_vector_view (x_init, p);

  gsl_rng * r;

  gsl_rng_env_setup();
  r = gsl_rng_alloc (gsl_rng_default);

  /* This is the data to be fitted */

  for (i = 0; i < n; i++)
    {
      double t = i;
      y[i] = 1.0 + 5 * exp (-0.1 * t) + gsl_ran_gaussian(r, 0.1);
      sigma[i] = 0.1;
      printf("data: %d %g %g\n", i, y[i], sigma[i]);
    };


  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, &f, &x);

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

      printf ("status = %s\n", gsl_strerror (status));

      print_state (iter, s);

      if (status)
	break;

      status = gsl_multifit_test_delta (s->dx, s->x, 0.0001, 0.0001);
    }
  while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);

  gsl_matrix_fprintf (stdout, covar, "%g");

#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  printf("b      = %g +/- %g\n", FIT(0), ERR(0));
  printf("A      = %g +/- %g\n", FIT(1), ERR(1));
  printf("lambda = %g +/- %g\n", FIT(2), ERR(2));


  printf ("status = %s\n", gsl_strerror (status));

  gsl_multifit_fdfsolver_free (s);
}

int
print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f |f(x)| = %g\n",
	  iter,
	  gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2), gsl_blas_dnrm2 (s->f));
}



