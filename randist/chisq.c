#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_chisq (const gsl_rng * r, double nu)
{
  double chisq = 2 * gsl_ran_gamma (r, nu / 2) ;
  return chisq ;
}

