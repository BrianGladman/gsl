#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_exponential (const gsl_rng * r, double mu)
{
  double u = gsl_rng_uniform_pos (r);

  return -mu * log (u);
}
