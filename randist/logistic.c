#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The logistic distribution has the form,

   p(x) dx = (1/mu) exp(-x/mu) / (1 + exp(-x/mu))^2 dx

   for x > 0 */

double
gsl_ran_logistic (const gsl_rng * r, const double mu)
{
  double x, z;

  do
    {
      x = gsl_rng_uniform_pos (r);
    }
  while (x == 1);

  z = log (x / (1 - x));

  return mu * z;
}

double
gsl_ran_logistic_pdf (const double x, const double mu)
{
  double u = exp (-x/mu);
  double p = u / (fabs(mu) * (1 + u) * (1 + u));
  return p;
}
