#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl_errno.h>
#include <gsl_rng.h>

/* gsl_rng_internal is a non-const version of gsl_rng, for creating gsl_rngs */

typedef struct {
  const char * name ;
  unsigned long int max ;
  size_t size ;
  void * state ;
  void (* set)(void * state, unsigned int seed) ;
  unsigned long int (* get)(void * state) ;
} gsl_rng_internal ;   

static gsl_rng * reference_to_generator (gsl_rng_internal * r) ;

gsl_rng *
gsl_rng_alloc (const gsl_rng_type * (*f)(void))
{
  
  const gsl_rng_type * T = f() ;

  gsl_rng_internal * r = (gsl_rng_internal *) malloc(sizeof(gsl_rng)) ;

  if (r == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for rng struct",
			GSL_ENOMEM, 0);
    } ;

  r->state = malloc(T->size) ;

  if (r->state == 0) 
    {
      free(r) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for rng state",
			GSL_ENOMEM, 0);
    } ;

  r->name = T->name ;
  r->max = T->max ;
  r->size = T->size ;
  r->set = T->set ;
  r->get = T->get ;

  gsl_rng_set ((gsl_rng *)r, gsl_rng_default_seed) ; /* seed the generator */

  return reference_to_generator (r);
}

static gsl_rng * reference_to_generator (gsl_rng_internal * r)
{
  gsl_rng * r_const = (gsl_rng *) r ;
  return r_const ;
}

void gsl_rng_set (const gsl_rng * r, unsigned int seed)
{
  (r->set)(r->state, seed) ;
}

unsigned long int gsl_rng_get (const gsl_rng * r)
{
  return (r->get)(r->state) ;
}

unsigned long int gsl_rng_max (const gsl_rng * r)
{
  return r->max ;
}

const char * gsl_rng_name (const gsl_rng * r)
{
  return r->name ;
}

void gsl_rng_print_state (const gsl_rng * r)
{
  size_t i ;
  unsigned char * p = (unsigned char *)(r->state) ; 
  const size_t n = r->size ;

  for (i = 0 ; i < n ; i++) 
    {
      printf("%.2x", *(p + i)) ; /* FIXME: assumed that a char is 8 bits */
    }
  
}

void gsl_rng_free (gsl_rng * r)
{
  free(r->state) ;
  free(r) ;
}


