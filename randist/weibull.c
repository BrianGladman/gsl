#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Weibull distribution has the form,

   p(x) dx = a x^(a-1) exp(-x^a) dx

 */

double
gsl_ran_weibull (const gsl_rng * r, const double mu, const double a)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (-log (x), 1 / a);

  return mu * z;
}

double
gsl_ran_weibull_pdf (const double x, const double mu, const double a)
{
  if (x < 0)
    {
      return 0 ;
    }
  else if (x == 0)
    {
      if (a == 1)
	return 1/mu ;
      else
	return 0 ;
    }
  else if (a == 1)
    {
      return exp(-x/mu)/mu ;
    }
  else
    {
      double p = (a/mu) * exp (-pow (x/mu, a) + (a - 1) * log (x/mu));
      return p;
    }
}
