#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Pareto distribution has the form,

   p(x) = a / x^(a+1)      for x >= 1

 */

double
gsl_ran_pareto (const gsl_rng * r, const double a)
{
  double x = gsl_rng_uniform_pos (r);

  double z = pow (x, -1 / a);

  return z;
}

double
gsl_ran_pareto_pdf (const double x, const double a)
{
  if (x >= 1)
    {
      double p = a / pow (x, a + 1);
      return p;
    }
  else
    {
      return 0;
    }
}

