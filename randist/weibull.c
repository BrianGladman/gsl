#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

/* The Weibull distribution is defined as,

   p(x) = b x^(b-1) exp(-x^b)

   */

double
gsl_ran_weibull (const gsl_rng * r, double a)
{
  double x = gsl_rng_uniform_pos (r);
  
  double z = pow(-log (x), 1/a) ;

  return z ;
}
