#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_diff.h>

double f (double x, void * params)
{
  return pow (x, 1.5);
}

int
main ()
{
  gsl_function F;
  double result, abserr;

  F.function = &f ;
  F.params = 0;

  printf("f(x) = x^(3/2)\n\n") ;

  gsl_diff_central (&F, 2.0, &result, &abserr);
  printf("x = 2.0\n");
  printf("derivative = %.10f +/- %.5f\n", result, abserr);
  printf("exact      = %.10f\n\n", 1.5 * sqrt(2.0));

  gsl_diff_forward (&F, 0.0, &result, &abserr);
  printf("x = 0.0\n");
  printf("derivative = %.10f +/- %.5f\n", result, abserr);
  printf("exact      = %.10f\n", 0.0);
}
