#include <config.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Geometric distribution (bernoulli trial with probability p) 

   prob(k) =  p (1 - p)^(k-1) for n = 1, 2, 3, ...

   It gives the distribution of "waiting times" for an event that
   occurs with probability p. */

unsigned int
gsl_ran_geometric (const gsl_rng * r, const double p)
{
  double u = gsl_rng_uniform_pos (r);

  unsigned int k;

  if (p == 1)
    {
      k = 1;
    }
  else
    {
      k = log (u) / log (1 - p) + 1;
    }

  return k;
}

double
gsl_ran_geometric_pdf (const unsigned int k, const double p)
{
  if (k == 0)
    {
      return 0 ;
    }
  else if (k == 1)
    {
      return p ;
    }
  else
    {
      double P = p * pow (1 - p, k - 1.0);
      return P;
    }
}
