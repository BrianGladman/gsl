#include <config.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

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
      double p, gab, ga, gb;
      gsl_sf_result result_ab, result_a, result_b ;

      gsl_sf_lngamma_impl (a + b, &result_ab);
      gsl_sf_lngamma_impl (a, &result_a);
      gsl_sf_lngamma_impl (b, &result_b);

      gab = result_ab.val ;
      ga = result_a.val ;
      gb = result_b.val ;
      
      p = exp (gab - ga - gb) * pow (x, a - 1) * pow (1 - x, b - 1);

      return p;
    }
}
