#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The logistic distribution has the form,

   p(x) dx = exp(-x) / (1 + exp(-x))^2 

   for x > 0 */

double
gsl_ran_logistic (const gsl_rng * r)
{
  double x, z;

  do
    {
      x = gsl_rng_uniform_pos (r);
    }
  while (x == 1);

  z = log (x / (1 - x));

  return z;
}

double
gsl_ran_logistic_pdf (double x)
{
  double u = exp (-x);
  double p = u / ((1 + u) * (1 + u));
  return p;
}
