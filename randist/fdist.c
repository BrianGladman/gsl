#include <config.h>
#include <math.h>
#include <gsl_sf_gamma.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The F distribution has the form

   p(x) dx = (nu1^(nu1/2) nu2^(nu2/2) Gamma((nu1 + nu2)/2) /
   Gamma(nu1/2) Gamma(nu2/2)) *
   x^(nu1/2 - 1) (nu2 + nu1 * x)^(-nu1/2 -nu2/2) dx

   The method used here is the one described in Knuth */

double
gsl_ran_fdist (const gsl_rng * r, const double nu1, const double nu2)
{

  double Y1 =  gsl_ran_gamma (r, nu1 / 2, 2.0);
  double Y2 =  gsl_ran_gamma (r, nu2 / 2, 2.0);

  double f = (Y1 * nu2) / (Y2 * nu1);

  return f;
}

double
gsl_ran_fdist_pdf (const double x, const double nu1, const double nu2)
{
  if (x < 0)
    {
      return 0 ;
    }
  else
    {
      double p;
      double lglg = (nu1 / 2) * log (nu1) + (nu2 / 2) * log (nu2) ;

      gsl_sf_result lg12, lg1, lg2 ;

      gsl_sf_lngamma_impl ((nu1 + nu2) / 2, &lg12);
      gsl_sf_lngamma_impl (nu1 / 2, &lg1);
      gsl_sf_lngamma_impl (nu2 / 2, &lg2);
      
      p = exp (lglg + lg12.val - lg1.val - lg2.val)
	* pow (x, nu1 / 2 - 1) * pow (nu2 + nu1 * x, -nu1 / 2 - nu2 / 2);
      
      return p;
    }
}
