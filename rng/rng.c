#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

gsl_rng *
gsl_rng_alloc (gsl_rng_type * T)
{
  gsl_rng * r ;
  r = (gsl_rng *) malloc(sizeof(gsl_rng)) ;
  r->state = malloc(T->state) ;
  r->set = T->set ;
  r->get = T->get ;
  gsl_rng_set(r, 1) ; /* seed the generator with a value of 1 */
  return r;
}

void gsl_rng_set (gsl_rng * r, unsigned int seed)
{
  (r->set)(r->state, seed) ;
}

unsigned long int gsl_rng_get (gsl_rng * r)
{
  return (r->get)(r->state) ;
}

void gsl_rng_free (gsl_rng * r)
{
  free(r->state) ;
  free(r) ;
}


