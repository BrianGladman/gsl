#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "demof.h"
#include "demof.c"

int
main ()
{
  int status;
  int iterations = 0, max_iterations = 100;
  gsl_root_fdfsolver *s;
  double x0, x = 5.0, r_expected = sqrt (5.0);
  gsl_function_fdf FDF;
  struct quadratic_params params = {1.0, 0.0, -5.0};

  FDF.f = &quadratic;
  FDF.df = &quadratic_deriv;
  FDF.fdf = &quadratic_fdf;
  FDF.params = &params;

  s = gsl_root_fdfsolver_alloc (gsl_root_fdfsolver_newton, &FDF, x);

  printf ("using %s method\n", gsl_root_fdfsolver_name (s));

  printf ("%-5s %10s %10s %10s %10s\n",
	  "iter", "root", "actual", "err", "err(est)");
  do
    {
      iterations++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, 0.001);

      if (status == GSL_SUCCESS)
	printf ("Converged:\n");

      printf ("%5d %10.7f %10.7f %+10.7f %10.7f\n",
	      iterations, x, r_expected, x - r_expected, x - x0);
    }
  while (status == GSL_CONTINUE && iterations < max_iterations);

}
