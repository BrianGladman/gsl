#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* Geometric distribution (bernoulli trial with probability p) 

   prob(n) =  p (1 - p)^(n-1) for n = 1, 2, 3, ...

   It gives the distribution of "waiting times" for an event that
   occurs with probability p. */

unsigned int
gsl_ran_geometric (const gsl_rng * r, const double p)
{
  double u = gsl_rng_uniform_pos (r);

  unsigned int n;

  if (p == 1)
    {
      n = 1;
    }
  else
    {
      n = log (u) / log (1 - p) + 1;
    }

  return n;
}

double
gsl_ran_geometric_pdf (const unsigned int n, const double p)
{
  if (n == 0)
    {
      return 0 ;
    }
  else if (n == 1)
    {
      return p ;
    }
  else
    {
      double P = p * pow (1 - p, n - 1.0);
      return P;
    }
}
