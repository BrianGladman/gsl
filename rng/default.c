#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl_rng.h>
#include <gsl_errno.h>

unsigned long int gsl_rng_default_seed = 1 ;
const gsl_rng_type * gsl_rng_default ;

static void 
check (const gsl_rng_type ** def, const gsl_rng_type * test, const char * p) ;

int
gsl_rng_env_setup (void)
{
  unsigned long int seed = 1;
  const char * p = getenv ("GSL_RNG_TYPE") ;

  if (p)
    {
      gsl_rng_default = 0 ;

      /* check GSL_RNG_TYPE against the names of all the generators */
      
      check(&gsl_rng_default, gsl_rng_bad_rand, p) ;
      check(&gsl_rng_default, gsl_rng_bad_randu, p) ;
      check(&gsl_rng_default, gsl_rng_cmrg, p) ;
      check(&gsl_rng_default, gsl_rng_mrg, p) ;
      check(&gsl_rng_default, gsl_rng_minstd, p) ;
      check(&gsl_rng_default, gsl_rng_rand, p) ;
      check(&gsl_rng_default, gsl_rng_taus, p) ;
      check(&gsl_rng_default, gsl_rng_uni, p) ;
      check(&gsl_rng_default, gsl_rng_uni32, p) ;
      check(&gsl_rng_default, gsl_rng_zuf, p) ;
      
      if (gsl_rng_default == 0) 
	{
	  GSL_ERROR("unknown generator",EINVAL) ;
	}
      
      printf("GSL_RNG_TYPE=%s\n",gsl_rng_default->name) ;
    }
  else
    {
      gsl_rng_default = gsl_rng_taus ;
    }

  p = getenv ("GSL_RNG_SEED") ;
   
   if (p) 
     {
       seed = strtoul(p, 0, 0) ;
       printf("GSL_RNG_SEED=%lu\n", seed) ;
     } ;

   gsl_rng_default_seed = seed ;

   return GSL_SUCCESS ;
}


static void 
check (const gsl_rng_type ** def, const gsl_rng_type * test, const char * p)
{
  if (*def)
    return ;

  if (strcmp(p,test->name) == 0) 
    {
      *def = test ;
    }
}


