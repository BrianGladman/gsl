#include <config.h>
#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* Logarithmic distribution 

   prob(n) =   p^n / (n log(1/(1-p)) for n = 1, 2, 3, ...

   We use Kemp's second accelerated generator, from Luc Devroye's book
   on "Non-Uniform Random Variate Generation", Springer */

unsigned int
gsl_ran_logarithmic (const gsl_rng * r, const double p)
{
  double c = log (1-p) ;

  double v = gsl_rng_uniform_pos (r);
  
  if (v >= p)
    {
      return 1 ;
    }
  else
    {
      double u = gsl_rng_uniform_pos (r);      
      double q = 1 - exp (c * u);

      if (v <= q*q)
	{
	  double x = 1 + log(v)/log(q) ;
	  return x ;
	}
      else if (v <= q)
	{
	  return 2;
	}
      else
	{
	  return 1 ;
	}
    }
}

double
gsl_ran_logarithmic_pdf (const unsigned int k, const double p)
{
  if (k == 0)
    {
      return 0 ;
    }
  else 
    {
      double P = pow(p, (double)k) / (double) k / log(1/(1-p)) ;
      return P;
    }
}
