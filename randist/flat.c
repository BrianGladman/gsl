#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* This is the uniform distribution in the range [a, b)

   p(x) dx = 1/(b-a) dx   if  a <= x < b
   .....   = 0            otherwise 

 */

double
gsl_ran_flat (const gsl_rng * r, const double a, const double b)
{
  double u = gsl_rng_uniform (r);

  /* A uniform distribution over [a,b] */

  return a * (1 - u) + b * u;
}

double
gsl_ran_flat_pdf (double x, const double a, const double b)
{
  if (x < b && x >= a)
    {
      return 1 / (b - a);
    }
  else
    {
      return 0;
    }
}
