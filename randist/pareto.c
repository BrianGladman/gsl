#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Pareto distribution has the form,

   p(x) dx = (a/mu) / (x/mu)^(a+1) dx     for x >= mu

 */

double
gsl_ran_pareto (const gsl_rng * r, double mu, const double a)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (x, -1 / a);

  return mu * z;
}

double
gsl_ran_pareto_pdf (const double x, const double mu, const double a)
{
  if (x >= mu)
    {
      double p = (a/mu) / pow (x/mu, a + 1);
      return p;
    }
  else
    {
      return 0;
    }
}

