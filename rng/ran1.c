#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is an implementation of the algorithm used in Numerical
   Recipe's ran1 generator.  It is MINSTD with a 32-element
   shuffle-box. 

   As far as I can tell, in general the effects of adding a shuffle
   box cannot be proven theoretically, so the period of this generator
   is unknown.  */

unsigned long int ran1_get (void *vstate);
void ran1_set (void *state, unsigned long int s);

static const long int m = 2147483647, a = 16807, q = 127773, r = 2836;

#define N_SHUFFLE 32
#define N_DIV (1 + 2147483646/N_SHUFFLE)

typedef struct
  {
    unsigned long int x;
    unsigned long int n;
    unsigned long int shuffle[N_SHUFFLE] ;
  }
ran1_state_t;

unsigned long int
ran1_get (void *vstate)
{
  ran1_state_t *state = (ran1_state_t *) vstate;

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

  { 
    unsigned long int j = state->n / N_DIV ;
    state->n = state->shuffle[j] ;
    state->shuffle[j] = state->x ;
  }

  return state->n ;
}


void
ran1_set (void *vstate, unsigned long int s)
{
  ran1_state_t *state = (ran1_state_t *) vstate;
  int i;

  if (s == 0)
    s = 1;	/* default seed is 1 */
  
  for (i = 0; i < 8; i++)
    {
      long int h = s / q;
      long int t = a * (s - h * q) - h * r;
      if (t < 0) t += m;
      s = t ;
    }

  for (i = N_SHUFFLE - 1; i >= 0; i--)
    {
      long int h = s / q;
      long int t = a * (s - h * q) - h * r;
      if (t < 0) t += m;
      s = t ;
      state->shuffle[i] = s ;
    }

  state->x = s ;
  state->n = s ;
 
  return;
}

static const gsl_rng_type ran1_type =
{"ran1",			/* name */
 2147483646,			/* RAND_MAX */
 1,         			/* RAND_MIN */
 sizeof (ran1_state_t),
 &ran1_set,
 &ran1_get};

const gsl_rng_type *gsl_rng_ran1 = &ran1_type;
