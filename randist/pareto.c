#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Pareto distribution is defined as,

   p(x) = a / x^(a+1)      for x >= 1

   */

double
gsl_ran_pareto (const gsl_rng * r, double a)
{
  double x = gsl_rng_uniform_pos (r);
  
  double z = pow (x, -1/a) ;

  return z ;
}
