#include <math.h>
#include <gsl_sf.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The sum of N samples from an exponential distribution gives an
   Erlang distribution

   p(x) = x^(n-1) exp (-x/a) / ((n-1)!a^n)

   for x = 0 ... +infty */

double
gsl_ran_erlang (const gsl_rng * r, double mu, double n)
{
  return mu * gsl_ran_gamma (r, n);
}

double
gsl_ran_erlang_pdf (double x, double mu, double n)
{
  if (x <= 0) 
    {
      return 0 ;
    }
  else
    {
      double lngamma = gsl_sf_lngamma (n);
      double p = exp ((n - 1) * log (x/mu) - x/mu - lngamma) / mu;
      return p;
    }
}
