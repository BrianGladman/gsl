#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

#include "seed.c"

/*  */

unsigned long int r250_get (void * vstate);
void r250_set (void * state, unsigned long int s);

typedef struct {
  int i ;
  unsigned long x[250];
} r250_state_t ;

unsigned long int r250_get (void *vstate)
{
    r250_state_t * state = (r250_state_t *)vstate;
    unsigned long int k ;
    int j ;

    int i = state->i ;

    if (i >= 147)
      {
	 j = i - 147;
      }
    else
      {
	j = i + 103;
      }
    
    k  = state->x[i] ^ state->x[j];
    state->x[i] = k;

    if (i >= 249 ) 
      {
	state->i = 0;
      }
    else
      {
	state->i = i + 1;
      }

    return k;
}

void r250_set(void * vstate, unsigned long int s)
{
  r250_state_t * state = (r250_state_t *) vstate;

  int i ;

  if (s == 0) s = 1 ; /* default seed is 1 */
  
  state->i = 0 ;

  for (i = 0; i < 250; i++)     /* Fill the buffer  */
    state->x[i] = lcg_seed(&s);

  for (i = 0; i < 250; i++)     /* Set some of the MS bits to 1 */
    if (lcg_seed(&s) > 0x20000000L)
      state->x[i] |= 0x40000000L;

  {
    /* Masks for turning on the diagonal bit and turning off the
       leftmost bits */

    unsigned long int msb =  0x40000000L; 
    unsigned long int mask = 0x7fffffffL; 
    
    for (i = 0; i < 31; i++)
      {
	int k = 7 * i + 3;         /* Select a word to operate on        */
	state->x[k] &= mask;       /* Turn off bits left of the diagonal */
	state->x[k] |= msb;        /* Turn on the diagonal bit           */
	mask >>= 1;
	msb >>= 1;
      }
  }

  return;
}

static const gsl_rng_type r250_type = { "r250",  /* name */
				        2147483647,  /* RAND_MAX */
					sizeof(r250_state_t), 
					&r250_set, 
					&r250_get } ;

const gsl_rng_type * gsl_rng_r250 = &r250_type ;
