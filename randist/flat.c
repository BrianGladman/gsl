#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_flat (const gsl_rng * r, double a, double b)
{
  double u = gsl_rng_uniform (r);

  /* A uniform distribution over [a,b] */

  return a * (1 - u) + b * u;
}
