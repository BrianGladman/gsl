#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is a reincarnation of the infamously bad RANDU generator,

   y -> (65539 * y) mod 2^31

   Knuth describes it as "really horrible".  */

unsigned long int bad_rand_get (void * vstate);
void bad_rand_set (void * state, unsigned int s);
void bad_rand_set_with_state (void * vstate, const void * vinit_state,
			 unsigned int s);

static const int a = 65539 ;
static const unsigned int m = 2147483648U ;
static const int q = 32766 ;
static const int r = 32774 ;

typedef struct {
  unsigned int x;
} bad_randu_state_t ;

unsigned long int bad_randu_get (void *vstate)
{
    bad_randu_state_t * state = (bad_randu_state_t *)vstate;
    int x = state->x ;

    unsigned int h = x / q ;
    x = q * (x - h * q) - h * r ;
    if (x < 0) x += m ;

    state->x = x ;

    return state->x;
}


void bad_randu_set(void * vstate, unsigned int s)
{
  bad_randu_state_t * state = (bad_randu_state_t *) vstate;
  
  state->x = s ;

  return;
}

static const gsl_rng_type bad_randu_type = { "bad_randu",  /* name */
					     2147483648UL,  /* RAND_MAX */
					     sizeof(bad_randu_state_t), 
					     &bad_randu_set, 
					     &bad_randu_get } ;

const gsl_rng_type * gsl_rng_bad_randu (void) { return &bad_randu_type ; }
