#include <config.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_sf.h>

/* The binomial distribution has the form,

   prob(k) =  n!/(k!(n-k)!) *  p^k (1-p)^(n-k) for k = 0, 1, ..., n

   This is the algorithm from Knuth */

unsigned int
gsl_ran_binomial (const gsl_rng * r, double p, unsigned int n)
{
  unsigned int i, a, b, k = 0;

  while (n > 10)	/* This parameter is tunable */
    {
      double X;
      a = 1 + (n / 2);
      b = 1 + n - a;

      X = gsl_ran_beta (r, (double) a, (double) b);

      if (X >= p)
	{
	  n = a - 1;
	  p /= X;
	}
      else
	{
	  k += a;
	  n = b - 1;
	  p = (p - X) / (1 - X);
	}
    }

  for (i = 0; i < n; i++)
    {
      double u = gsl_rng_uniform (r);
      if (u < p)
	k++;
    }

  return k;
}

double
gsl_ran_binomial_pdf (const unsigned int k, const double p, 
		      const unsigned int n)
{
  if (k > n)
    {
      return 0 ;
    }
  else 
    {
      double a = k;
      double b = n - k;
      double P;
      gsl_sf_result cnk ;

      gsl_sf_choose_impl (n, k, &cnk) ;

      P = cnk.val * pow (p, a) * pow (1 - p, b);
      
      return P;
    }
}
