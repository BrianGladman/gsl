#include <config.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The exponential distribution has the form

   p(x) dx = exp(-x/mu) dx/mu

   for x = 0 ... +infty */

double
gsl_ran_exponential (const gsl_rng * r, const double mu)
{
  double u = gsl_rng_uniform_pos (r);

  return -mu * log (u);
}

double
gsl_ran_exponential_pdf (const double x, const double mu)
{
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      double p = exp (-x/mu)/mu;
      
      return p;
    }
}
