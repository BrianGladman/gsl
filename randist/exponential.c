#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_exponential (const gsl_rng * r, double mu)
{
  double u;

  do
    {
      u = gsl_rng_uniform (r);
    }
  while (u == 0.0);

  return -mu * log (u);
}
