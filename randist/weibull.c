#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Weibull distribution has the form,

   p(x) = a x^(a-1) exp(-x^a)

 */

double
gsl_ran_weibull (const gsl_rng * r, const double a)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (-log (x), 1 / a);

  return z;
}

double
gsl_ran_weibull_pdf (const double x, const double a)
{
  if (x < 0)
    {
      return 0 ;
    }
  else if (x == 0)
    {
      if (a == 1)
	return 1 ;
      else
	return 0 ;
    }
  else if (a == 1)
    {
      return exp(-x) ;
    }
  else
    {
      double p = a * exp (-pow (x, a) + (a - 1) * log (x));
      return p;
    }
}
