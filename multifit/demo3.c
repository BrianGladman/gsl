#include <stdio.h>
#include <math.h>
#include <gsl/gsl_randist.h>

int main (void)
{
  double x;
  const gsl_rng_type * T;
  gsl_rng * r;
  
  gsl_rng_env_setup();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  for (x = 0.1 ; x < 2; x+= 0.1)
    {
      double y0 = exp(x);
      double err = 0.1*y0;

      printf("%g %g %g\n", x, y0 + gsl_ran_gaussian(r, err), err);
    }
}
  
