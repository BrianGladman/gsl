#include <stdio.h>
#include <gsl_errno.h>
#include <gsl_math.h>
#include <gsl_min.h>

double fn1 (double x, void * params);

double fn1 (double x, void * params)
{
  return cos(x) + 1.0 ;
}

int
main ()
{
  int status;
  int iterations = 0, max_iterations = 100;
  gsl_min_fminimizer *s;
  double m = 2.0, m_expected = M_PI;
  gsl_interval x =  {0, 6.0};
  gsl_function F;

  F.function = &fn1;
  F.params = 0;

  s = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent, &F, m, x);

  printf ("using %s method\n", gsl_min_fminimizer_name (s));

  printf ("%5s [%9s, %9s] %9s %9s %10s %9s\n",
	  "iter", "lower", "upper", "min", "actual", "err", "err(est)");

  printf ("%5d [%.7f, %.7f] %.7f %.7f %+.7f %.7f\n",
          iterations, x.lower, x.upper,
          m, m_expected, m - m_expected, x.upper - x.lower);

  do
    {
      iterations++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_minimum (s);
      x = gsl_min_fminimizer_interval (s);

      status = gsl_min_test_interval (x, 0.001, 0.0);

      if (status == GSL_SUCCESS)
	printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %.7f %+.7f %.7f\n",
	      iterations, x.lower, x.upper,
	      m, m_expected, m - m_expected, x.upper - x.lower);
    }
  while (status == GSL_CONTINUE && iterations < max_iterations);

}
