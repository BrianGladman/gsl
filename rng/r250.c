#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

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
  unsigned long int * const x = state->x ;

  state->i = 0 ;

  for (i = 0; i < 250; i++)     /* Fill the buffer with 15-bit values */
    state->x[j] = rand();

  for (j = 0; j < 250; j++)     /* Set some of the MS bits to 1 */
    if (rand() > 0x20000000L)
      r250_buffer[j] |= 0x40000000L;

  msb =  0x40000000L;      /* To turn on the diagonal bit   */
  mask = 0x7fffffffL;      /* To turn off the leftmost bits */

  for (j = 0; j < 31; j++)
    {
      int k = 7 * j + 3;         /* Select a word to operate on        */
      x[k] &= mask;              /* Turn off bits left of the diagonal */
      x[k] |= msb;               /* Turn on the diagonal bit           */
      mask >>= 1;
      msb >>= 1;
    }
      
  return;
}

static const gsl_rng_type r250_type = { "r250",  /* name */
				       4294967295UL,  /* RAND_MAX */
				       sizeof(r250_state_t), 
				       &r250_set, 
				       &r250_get } ;

const gsl_rng_type * gsl_rng_r250 = &r250_type ;
