#include <math.h>
#include <gsl_sf.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The chisq distributions has the form

   p(x) dx = (1/Gamma(nu/2)) (x/2)^(nu/2 - 1) exp(-x/2) dx

   for x = 0 ... +infty */

double
gsl_ran_chisq (const gsl_rng * r, const double nu)
{
  double chisq = 2 * gsl_ran_gamma (r, nu / 2, 1.0);
  return chisq;
}

double
gsl_ran_chisq_pdf (const double x, const double nu)
{
  if (x <= 0)
    {
      return 0 ;
    }
  else
    {
      double lngamma = gsl_sf_lngamma (nu / 2);
      double p = exp ((nu / 2 - 1) * log (x/2) - x/2 - lngamma) / 2;
      return p;
    }
}
