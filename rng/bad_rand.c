#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is the example generator given in the ANSI C standard. It is
   not very good. */

unsigned long int bad_rand_get (void * vstate);
void bad_rand_set (void * state, unsigned int s);
void bad_rand_set_with_state (void * vstate, const void * vinit_state,
			 unsigned int s);

static const int a = 1103515245 ;
static const int c = 12345 ;
static const int m = 32768 ;

typedef struct {
  long int x;
} bad_rand_state_t ;

unsigned long int bad_rand_get (void *vstate)
{
    bad_rand_state_t * state = (bad_rand_state_t *)vstate;

    unsigned long int x = (a * state->x) + c ;
    
    state->x = (x / 65536) % m ;

    return state->x;
}


void bad_rand_set(void * vstate, unsigned int s)
{
  bad_rand_state_t * state = (bad_rand_state_t *) vstate;
  
  state->x = s ;

  return;
}

static const gsl_rng_type bad_rand_type = { "bad_rand",  /* name */
					     32768,  /* RAND_MAX */
					     sizeof(bad_rand_state_t), 
					     &bad_rand_set, 
					     &bad_rand_get } ;

const gsl_rng_type * gsl_rng_bad_rand (void) { return &bad_rand_type ; }
