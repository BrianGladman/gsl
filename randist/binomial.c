#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* Binomial distribution 

   prob(t) =  (t!/(n!(t-n)!)) *  p^n (1 - p)^(N-n) for n = 0, 1, ..., N

   */

unsigned int
gsl_ran_binomial (const gsl_rng * r, double p, unsigned int t)
{
  p = 0;
  r = 0;
  t = 0;

  return 0;
}
