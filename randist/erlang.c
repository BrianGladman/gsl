#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The sum of N samples from an exponential distribution gives an
   Erlang distribution 

   p(x) = x^(n-1) exp (-x/a) / ((n-1)!a^n)

   */

double
gsl_ran_erlang (const gsl_rng * r, double mu, unsigned int n)
{
  return mu * gsl_ran_gamma(r, (double)n);
}
