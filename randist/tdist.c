#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_tdist (const gsl_rng * r, double nu)
{
  /* The t-distribution has the form

     p(x) dx = (Gamma((nu + 1)/2)/(sqrt(pi nu) Gamma(nu/2)) *
                 (1 + (x^2)/nu)^-((nu + 1)/2)

     The method used here is the one described in Knuth */

  if (nu <= 2)
    {
      double y1 = gsl_ran_gaussian (r) ;
      double y2 = gsl_ran_chisq (r, nu) ;
      
      double t = y1 / sqrt(y2 / nu) ;

      return t ;
    }
  else
    {
      double y1, y2, z, t ;
      do 
	{
	  y1 = gsl_ran_gaussian (r) ;
	  y2 = gsl_ran_exponential (r, 2/(nu - 2)) ;
	  
	  z = y1*y1 / (nu - 2) ;
	}
      while (exp(-y2-z) >= (1 - z)) ;
      
      t = y1 / sqrt((1-2*nu)*(1-z)) ;
      return t ;
    }
}
