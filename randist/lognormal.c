#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

#include "gaussian.c"

double
gsl_ran_lognormal (const gsl_rng * r)
{
  double x = gaussian (r);  /* First compute a gaussian */

  /* Lognormal numbers are the exponentials of gaussian random numbers */

  double z = exp (x) ;

  return z ;
}
