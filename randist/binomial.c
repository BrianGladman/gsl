#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_sf.h>

/* Binomial distribution 

   prob(t) =  (t!/(n!(t-n)!)) *  p^n (1 - p)^(t-n) for n = 0, 1, ..., t

   This is the algorithm from Knuth */

unsigned int
gsl_ran_binomial (const gsl_rng * r, double p, unsigned int t)
{
  unsigned int i, a, b, n = 0;

  while (t > 10)  /* This paramter is tunable */
    {
      double X ;
      a = 1 + (t / 2) ;
      b = 1 + t - a ;
      
      X = gsl_ran_beta (r, (double) a, (double) b) ;
      
      if (X >= p)
	{
	  t = a - 1 ;
	  p /= X ;
	}
      else
	{
	  n += a ;
	  t = b - 1 ;
	  p = (p - X)/(1 - X) ;
	}
    }

  for (i = 0; i < t ; i++)
    {
      double u = gsl_rng_uniform (r) ;
      if (u < p)
	n++ ;
    }

  return n;
}

double
gsl_ran_binomial_pdf (unsigned int n, double p, unsigned int t)
{
  double P = (gsl_sf_choose (t, n) 
	      * pow(p, (double)n) 
	      * pow(1 - p, (double)t - (double)n));
  return P ;
}
