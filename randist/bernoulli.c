#include <config.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

/* The bernoulli distribution has the form,

   prob(0) = 1-p, prob(1) = p

   */

unsigned int
gsl_ran_bernoulli (const gsl_rng * r, double p)
{
  double u = gsl_rng_uniform (r) ;

  if (u < p)
    {
      return 1 ;
    }
  else
    {
      return 0 ;
    }
}

double
gsl_ran_bernoulli_pdf (const unsigned int k, double p)
{
  if (k == 0)
    {
      return 1 - p ;
    }
  else if (k == 1)
    {
      return p ;
    }
  else
    {
      return 0 ;
    }
}
