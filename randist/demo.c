#include <stdio.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

int
main ()
{
  gsl_rng * r ;

  int i, n = 10;
  double mu = 3.0;

  /* create a generator chosen by the environment variable GSL_RNG_TYPE */

  gsl_rng_env_setup();
  
  r = gsl_rng_alloc (gsl_rng_default);

  /* print n random variates chosen from the poisson distribution with
     mean parameter mu */

  for (i = 0; i < n; i++) 
    {
      unsigned int k = gsl_ran_poisson (r, mu);
      
      printf(" %u", k);
    }

  printf("\n");
}
