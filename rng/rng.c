#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_rng.h>

gsl_rng *
gsl_rng_alloc (const gsl_rng_type * (*f)(void))
{
  gsl_rng * r ;
  const gsl_rng_type * T = f() ;

  r = (gsl_rng *) malloc(sizeof(gsl_rng)) ;

  if (r == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for rng struct",
			GSL_ENOMEM, 0);
    } ;

  r->state = malloc(T->state) ;


  if (r->state == 0) 
    {
      free(r) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for rng state",
			GSL_ENOMEM, 0);
    } ;


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


