#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_fdist (const gsl_rng * r, double nu1, double nu2)
{
  /* The F distribution has the form

     p(x) dx = (nu1^(nu1/2) nu2^(nu2/2) Gamma((nu1 + nu2)/2) /
                 Gamma(nu1/2) Gamma(nu2/2)) *
		  x^(nu1/2 - 1) (nu2 + nu1 * x)^(-nu1/2 -nu2/2)

     The method used here is the one described in Knuth */

  double Y1 = 2 * gsl_ran_gamma (r, nu1 / 2) ;
  double Y2 = 2 * gsl_ran_gamma (r, nu2 / 2) ;

  double f = (Y1 * nu2) / (Y2 * nu1) ;
  
  return f ;
}
