#include <stdio.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_roots.h>

#include "demof.h"
#include "demof.c"

int
main ()
{
  int status;
  int iterations = 0, max_iterations = 100;
  gsl_root_fsolver *s;
  double r = 0, r_expected = sqrt (5.0);
  gsl_interval x =
  {0.0, 5.0};
  gsl_function F;
  struct quadratic_params params =
  {1.0, 0.0, -5.0};

  F.function = &quadratic;
  F.params = &params;

  s = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection, &F, x);

  printf ("using %s method\n", gsl_root_fsolver_name (s));

  printf ("%5s [%9s, %9s] %9s %9s %10s %9s\n",
	  "iter", "lower", "upper", "root", "actual", "err", "err(est)");

  do
    {
      iterations++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      x = gsl_root_fsolver_interval (s);
      status = gsl_root_test_interval (x, 0, 0.001);

      if (status == GSL_SUCCESS)
	printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %.7f %+.7f %.7f\n",
	      iterations, x.lower, x.upper,
	      r, r_expected, r - r_expected, x.upper - x.lower);
    }
  while (status == GSL_CONTINUE && iterations < max_iterations);

}
