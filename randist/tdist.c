#include <math.h>
#include <gsl_math.h>
#include <gsl_sf.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The t-distribution has the form

   p(x) dx = (Gamma((nu + 1)/2)/(sqrt(pi nu) Gamma(nu/2))
   * (1 + (x^2)/nu)^-((nu + 1)/2)

   The method used here is the one described in Knuth */

double
gsl_ran_tdist (const gsl_rng * r, double nu)
{
  if (nu <= 2)
    {
      double Y1 = gsl_ran_gaussian (r);
      double Y2 = gsl_ran_chisq (r, nu);

      double t = Y1 / sqrt (Y2 / nu);

      return t;
    }
  else
    {
      double Y1, Y2, Z, t;
      do
	{
	  Y1 = gsl_ran_gaussian (r);
	  Y2 = gsl_ran_exponential (r, 2 / (nu - 2));

	  Z = Y1 * Y1 / (nu - 2);
	}
      while (exp (-Y2 - Z) >= (1 - Z));

      /* FIXME: there must be a typo in Knuth's formula
	 sqrt(1-2nu) can't be right if nu > 2 */
      t = Y1 / sqrt ((1 - 2 * nu) * (1 - Z));
      return t;
    }
}

double
gsl_ran_tdist_pdf (double x, double nu)
{
  double lg2 = gsl_sf_lngamma ((nu + 1) / 2);
  double lg1 = gsl_sf_lngamma (nu / 2);

  double p = exp (lg2 - lg1) / sqrt (M_PI * nu) * pow ((1 + x * x / nu),
						       -(nu + 1) / 2);
  return p;
}
