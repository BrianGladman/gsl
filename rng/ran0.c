#include <config.h>
#include <stdlib.h>
#include <gsl_errno.h>
#include <gsl_rng.h>

/* This is an implementation of the algorithm used in Numerical
   Recipe's ran0 generator. It is the same as MINSTD with an XOR mask
   of 123459876.

   The period of this generator is 2^31.  

   Note, if you choose a seed of 123459876 it would give a degenerate
   series 0,0,0,0, ...  I've made that into an error. */

unsigned long int ran0_get (void *vstate);
void ran0_set (void *state, unsigned long int s);

static const long int m = 2147483647, a = 16807, q = 127773, r = 2836;
static const unsigned long int mask = 123459876 ;

typedef struct
  {
    unsigned long int x;
  }
ran0_state_t;

unsigned long int
ran0_get (void *vstate)
{
  ran0_state_t *state = (ran0_state_t *) vstate;

  const unsigned long int x = state->x;

  const long int h = x / q;
  const long int t = a * (x - h * q) - h * r;

  if (t < 0)
    {
      state->x = t + m;
    }
  else
    {
      state->x = t;
    }

  return state->x ;
}


void
ran0_set (void *vstate, unsigned long int s)
{
  ran0_state_t *state = (ran0_state_t *) vstate;

  if (s == mask)
    {
      GSL_ERROR_RETURN_NOTHING ("ran0 should not use seed == mask", GSL_EINVAL) ;
    }

  state->x = s ^ mask;

  return;
}

static const gsl_rng_type ran0_type =
{"ran0",			/* name */
 2147483646,			/* RAND_MAX */
 1,         			/* RAND_MIN */
 sizeof (ran0_state_t),
 &ran0_set,
 &ran0_get};

const gsl_rng_type *gsl_rng_ran0 = &ran0_type;
