#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Laplace probability distribution is  

   p(x) = (1/(2 mu)) * exp( -|x/mu|)

   for -infty < x < infty  */

double
gsl_ran_two_sided_exponential (const gsl_rng * r, double mu)
{
  double u;
  do
    {
      u = 2 * gsl_rng_uniform (r) - 1.0;
    }
  while (u == 0.0);

  if (u < 0)
    {
      return mu * log (-u);
    }
  else
    {
      return -mu * log (u);
    }
}
