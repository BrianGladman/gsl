#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_sf.h>

/* The bernoulli distribution has the form,

   prob(0) = 1-p, prob(1) = p

   */

unsigned int
gsl_ran_negative_binomial (const gsl_rng * r, double p, unsigned int t)
{
  double X = gsl_ran_gamma (r, t) ;
  unsigned int n = gsl_ran_poisson (r, X*(1-p)/p) ;
  return n ;
}

double
gsl_ran_negative_binomial_pdf (const unsigned int n, const double p, 
			       const unsigned int t)
{
  double T = t ;
  double N = n ;
  double P = gsl_sf_choose (T-1+N, N) * pow (p, T) * pow (1 - p, N);
  
  return P;
}
