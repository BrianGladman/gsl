#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Weibull distribution has the form,

   p(x) = a x^(a-1) exp(-x^a)

 */

double
gsl_ran_weibull (const gsl_rng * r, double a)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (-log (x), 1 / a);

  return z;
}

double
gsl_ran_weibull_pdf (double x, double a)
{
  double p = a * exp (-pow (x, a) + (a - 1) * log (x));
  return p;
}
