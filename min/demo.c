#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

double fn1 (double x, void * params);

double fn1 (double x, void * params)
{
  return cos(x) + 1.0 ;
}

int
main ()
{
  int status;
  int iter = 0, max_iter = 100;
  gsl_min_fminimizer *s;
  double m = 2.0, m_expected = M_PI;
  double x_lower =  0, x_upper = 6.0;
  gsl_function F;

  F.function = &fn1;
  F.params = 0;

  s = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
  gsl_min_fminimizer_set (s, &F, m, x_lower, x_upper);

  printf ("using %s method\n", gsl_min_fminimizer_name (s));

  printf ("%5s [%9s, %9s] %9s %9s %10s %9s\n",
          "iter", "lower", "upper", "min", "actual",
          "err", "err(est)");

  printf ("%5d [%.7f, %.7f] %.7f %.7f %+.7f %.7f\n",
          iter, x_lower, x_upper,
          m, m_expected, m - m_expected, x_upper - x_lower);

  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);

      m = gsl_min_fminimizer_x_minimum (s);
      x_lower = gsl_min_fminimizer_x_lower (s);
      x_upper = gsl_min_fminimizer_x_upper (s);

      status = gsl_min_test_interval (x_lower, x_upper, 0.001, 0.0);

      if (status == GSL_SUCCESS)
        printf ("Converged:\n");

      printf ("%5d [%.7f, %.7f] %.7f %.7f %+.7f %.7f\n",
              iter, x_lower, x_upper,
              m, m_expected, m - m_expected, x_upper - x_lower);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

}
