#include <math.h>
#include <gsl_sf.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The sum of N samples from an exponential distribution gives an
   Erlang distribution

   p(x) dx = x^(n-1) exp (-x/a) / ((n-1)!a^n) dx

   for x = 0 ... +infty */

double
gsl_ran_erlang (const gsl_rng * r, const double a, const double n)
{
  return a * gsl_ran_gamma (r, n);
}

double
gsl_ran_erlang_pdf (const double x, const double a, const double n)
{
  if (x <= 0) 
    {
      return 0 ;
    }
  else
    {
      double lngamma = gsl_sf_lngamma (n);
      double p = exp ((n - 1) * log (x/a) - x/a - lngamma) / a;
      return p;
    }
}
