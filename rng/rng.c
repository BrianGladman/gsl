#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl_errno.h>
#include <gsl_rng.h>

gsl_rng *
gsl_rng_alloc (const gsl_rng_type * T)
{
  
  gsl_rng * r = (gsl_rng *) malloc(sizeof(gsl_rng)) ;

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

  gsl_rng_set (r, gsl_rng_default_seed) ; /* seed the generator */

  return r;
}


gsl_rng *
gsl_rng_clone (const gsl_rng * q)
{
  gsl_rng * r = (gsl_rng *) malloc(sizeof(gsl_rng)) ;

  if (r == 0) 
    {
      GSL_ERROR_RETURN ("failed to allocate space for rng struct",
			GSL_ENOMEM, 0);
    } ;

  r->state = malloc(q->size) ;

  if (r->state == 0) 
    {
      free(r) ; /* exception in constructor, avoid memory leak */

      GSL_ERROR_RETURN ("failed to allocate space for rng state",
			GSL_ENOMEM, 0);
    } ;

  r->name = q->name ;
  r->max = q->max ;
  r->size = q->size ;
  r->set = q->set ;
  r->get = q->get ;

  memcpy(r->state, q->state, q->size) ;

  return r;
}

void gsl_rng_set (const gsl_rng * r, unsigned long int seed)
{
  (r->set)(r->state, seed) ;
}

unsigned long int gsl_rng_get (const gsl_rng * r)
{
  return (r->get)(r->state) ;
}

double gsl_rng_get_uni (const gsl_rng * r)
{
  unsigned long int k = (r->get)(r->state) ;
  unsigned long int max = r->max ;

  return k / (1.0 + max) ;
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


