#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl_rng.h>

void get_default_seed (void) ;

unsigned int gsl_rng_default_seed = 1 ;

void
get_default_seed (void)
{
  unsigned int seed = 1;

  const char * p = getenv ("GSL_RNG_SEED") ;
   
   if (p)  /* GSL_IEEE_MODE environment variable is not set */
     {
       seed = strtoul(p, 0, 0) ;
       printf("GSL_RNG_SEED=%d\n", seed) ;
     } ;

   gsl_rng_default_seed = seed ;
}
