#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>

double
f (double x, void *p)
{
  p = 0;
  if (x < 0.5)
    return 0.25;
  else
    return 0.75;
}


int
main (void)
{
  double x;

  gsl_cheb_series *cs = gsl_cheb_alloc (40);

  gsl_function F;

  F.function = f;
  F.params = 0;

  gsl_cheb_init (cs, &F, 0.0, 1.0);

  for (x = -1; x < 2; x += 0.001)
    {
      double r10 = gsl_cheb_eval_n (cs, 10, x);
      double r40 = gsl_cheb_eval (cs, x);
      printf ("%g %g %g %g\n", x, GSL_FN_EVAL (&F, x), r5, r40);
    }

  gsl_cheb_free (cs);
}
