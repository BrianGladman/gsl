#ifndef GSL_RNG_H
#define GSL_RNG_H
#include <stdlib.h>

typedef struct {
  size_t n ;
  unsigned long int max ;
  void * state ;
  void (* set)(void * state, unsigned int seed) ;
  unsigned long int (* get)(void * state) ;
} gsl_rng ;

typedef struct { 
  size_t state ;
  void (* set)(void * state, unsigned int seed) ;
  unsigned long int (* get)(void * state) ;
} gsl_rng_type ;

const gsl_rng_type * gsl_rng_blah (void) ;

gsl_rng * gsl_rng_alloc (const gsl_rng_type * (*f)(void)) ;
unsigned long int gsl_rng_get (gsl_rng * r) ;
void gsl_rng_set (gsl_rng * r, unsigned int seed) ;
void gsl_rng_free (gsl_rng * r) ;

#endif /* GSL_RNG_H */

