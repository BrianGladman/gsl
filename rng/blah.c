#include <config.h>
#include <gsl_rng.h>

void blah_set (void * state, unsigned int seed) ;
unsigned long int blah_get (void * state) ;

static const 
gsl_rng_type blah_type = { sizeof(double), 
			   &blah_set,
			   &blah_get } ;

const gsl_rng_type * 
gsl_rng_blah (void) 
{
  return &blah_type ;
}

void 
blah_set (void * state, unsigned int seed)
{
  *(double *)state = seed ;
}

unsigned long int 
blah_get (void * state)
{
  *(double *)state = *(double *)state + 1 ;
  return (unsigned long int)(*(double *)state) ;
}
