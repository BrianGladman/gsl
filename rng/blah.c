#include <config.h>
#include <gsl_rng.h>

void 
gsl_blah_set (double seed, void * state)
{
  *(double *)state = seed ;
}

unsigned long int 
gsl_blah_get (void * state)
{
  *(double *)state = *(double *)state + 1 ;
  return (unsigned long int)(*(double *)state) ;
}
