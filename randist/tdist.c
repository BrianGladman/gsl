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
      double Y1 = gsl_ran_gaussian (r) ;
      double Y2 = gsl_ran_chisq (r, nu) ;
      
      double t = Y1 / sqrt(Y2 / nu) ;

      return t ;
    }
  else
    {
      double Y1, Y2, Z, t ;
      do 
	{
	  Y1 = gsl_ran_gaussian (r) ;
	  Y2 = gsl_ran_exponential (r, 2/(nu - 2)) ;
	  
	  Z = Y1*Y1 / (nu - 2) ;
	}
      while (exp(-Y2-Z) >= (1 - Z)) ;
      
      t = Y1 / sqrt((1-2*nu)*(1-Z)) ;
      return t ;
    }
}
