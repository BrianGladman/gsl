#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <gsl_sf.h>

/* The beta distribution has the form

   p(x) dx = (Gamma(a + b)/(Gamma(a) Gamma(b))) x^(a-1) (1-x)^(b-1) dx

   The method used here is the one described in Knuth */

double
gsl_ran_beta (const gsl_rng * r, const double a, const double b)
{
  double x1 = gsl_ran_gamma (r, a, 1.0);
  double x2 = gsl_ran_gamma (r, b, 1.0);

  return x1 / (x1 + x2);
}

double
gsl_ran_beta_pdf (const double x, const double a, const double b)
{
  if (x < 0 || x > 1)
    {
      return 0 ;
    }
  else 
    {
      double gab = gsl_sf_lngamma (a + b);
      double ga = gsl_sf_lngamma (a);
      double gb = gsl_sf_lngamma (b);
      
      double p = exp (gab - ga - gb) * pow (x, a - 1) * pow (1 - x, b - 1);
      return p;
    }
}
