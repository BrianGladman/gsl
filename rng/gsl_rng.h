#ifndef GSL_RNG_H
#define GSL_RNG_H
#include <stdlib.h>

typedef struct {
  const char * const name ;
  const unsigned long int max ;
  size_t size ;
  void * const state ;
  void (* const set)(void * state, unsigned int seed) ;
  unsigned long int (* const get)(void * state) ;
} gsl_rng ;

typedef struct { 
  const char * name ;
  const unsigned long int max ;
  size_t size ;
  void (* set)(void * state, unsigned int seed) ;
  unsigned long int (* get)(void * state) ;
} gsl_rng_type ;

const gsl_rng_type * gsl_rng_blah (void) ;
const gsl_rng_type * gsl_rng_cmrg (void) ;
const gsl_rng_type * gsl_rng_mrg (void) ;
const gsl_rng_type * gsl_rng_rand (void) ;

extern unsigned int gsl_rng_default_seed ;

gsl_rng * gsl_rng_alloc (const gsl_rng_type * (*f)(void)) ;
unsigned long int gsl_rng_get (const gsl_rng * r) ;
void gsl_rng_set (const gsl_rng * r, unsigned int seed) ;
unsigned long int gsl_rng_max (const gsl_rng * r) ;
const char * gsl_rng_name (const gsl_rng * r) ;
void gsl_rng_print_state (const gsl_rng * r) ;
void gsl_rng_free (gsl_rng * r) ;

#endif /* GSL_RNG_H */

