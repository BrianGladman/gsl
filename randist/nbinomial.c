#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_sf.h>

/* The negative binomial distribution has the form,

   prob(k) =  Gamma(n + k)/(Gamma(n) Gamma(k + 1))  p^n (1-p)^k 

   for k = 0, 1, ... . Note that n does not have to be an integer.

   This is the Leger's algorithm (given in the answers in Knuth) */

unsigned int
gsl_ran_negative_binomial (const gsl_rng * r, double p, double n)
{
  double X = gsl_ran_gamma (r, n, 1.0) ;
  unsigned int k = gsl_ran_poisson (r, X*(1-p)/p) ;
  return k ;
}

double
gsl_ran_negative_binomial_pdf (const unsigned int k, const double p, double n)
{
  double f = gsl_sf_lngamma (k + n) ;
  double a = gsl_sf_lngamma (n) ;
  double b = gsl_sf_lngamma (k + 1.0) ;
  double P = exp(f-a-b) * pow (p, n) * pow (1 - p, (double)k);
  
  return P;
}
