#include <math.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

double
gsl_ran_logistic (const gsl_rng * r)
{
  double x, z ;

  do 
    {
      x = gsl_rng_uniform_pos (r);  
    }
  while (x == 1) ;

  z = log(x/(1-x)) ;

  return z ;
}
