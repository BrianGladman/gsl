#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_sf.h>

/* The Pascal distribution is a negative binomial with valued integer n

   prob(k) =  (n - 1 + k)!/(n!(k - 1)!) *  p^n (1-p)^k for k = 0, 1, ..., n

   */

unsigned int
gsl_ran_pascal (const gsl_rng * r, double p, unsigned int n)
{
  /* This is a separate interface for the pascal distribution so that
     it can be optimized differently from the negative binomial in
     future.
     
     e.g. if n < 10 it might be faster to generate the Pascal
     distributions as the sum of geometric variates directly.  */
  
  unsigned int k = gsl_ran_negative_binomial (r, p, (double) n);
  return k;
}

double
gsl_ran_pascal_pdf (const unsigned int k, const double p, unsigned int n)
{
  double P = gsl_ran_negative_binomial_pdf (k, p, (double) n);
  return P;
}
