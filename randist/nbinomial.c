#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_sf.h>

/* The negative binomial distribution has the form,

   prob(t) =  (t - 1 + n)!/(n!(t - 1)!) *  p^t (1-p)^n for n = 0, 1, ..., t

   This is the Leger's algorithm (given in the answers in Knuth) */

unsigned int
gsl_ran_negative_binomial (const gsl_rng * r, double p, double t)
{
  double X = gsl_ran_gamma (r, t) ;
  unsigned int n = gsl_ran_poisson (r, X*(1-p)/p) ;
  return n ;
}

double
gsl_ran_negative_binomial_pdf (const unsigned int n, const double p, double t)
{
  double f = gsl_sf_lngamma (t + n) ;
  double a = gsl_sf_lngamma (t) ;
  double b = gsl_sf_lngamma (n + 1.0) ;
  double P = exp(f-a-b) * pow (p, t) * pow (1 - p, (double)n);
  
  return P;
}
