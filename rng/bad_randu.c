#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is a reincarnation of the infamously bad RANDU generator.
   The sequence is,
   
         x_{n+1} = (a x_n) mod m

   with a = 65539 and m = 2^31 = 2147483648. The seed specifies
   the initial value, x_1.

   The theoretical value of x_{10001} is 1623524161.

   The period of this generator is 2^29.

   Note: Knuth describes this generator as "really horrible". 

   From: Park and Miller, "Random Number Generators: Good ones are
   hard to find" Communications of the ACM, October 1988, Volume 31,
   No 10, pages 1192-1201. */

unsigned long int bad_randu_get (void * vstate);
void bad_randu_set (void * state, unsigned long int s);

static const long int a = 65539 ;
static const unsigned long int m = 2147483648UL ;
static const long int q = 32766 ;
static const long int r = 32774 ;

typedef struct {
  unsigned long int x;
} bad_randu_state_t ;

unsigned long int 
bad_randu_get (void *vstate)
{
    bad_randu_state_t * state = (bad_randu_state_t *)vstate;

    const unsigned long int x = state->x ;

    const long int h = x / q ;
    const long int t = a * (x - h * q) - h * r ;

    if (t < 0) 
      {
	state->x = t + m;
      }
    else
      {
	state->x = t ;
      }
    
    return state->x;
}

void 
bad_randu_set(void * vstate, unsigned long int s)
{
  bad_randu_state_t * state = (bad_randu_state_t *) vstate;

  if (s == 0) s = 1 ; /* default seed is 1 */
  
  state->x = s ;

  return;
}

static const gsl_rng_type bad_randu_type = { "bad-randu",  /* name */
					     2147483648UL,  /* RAND_MAX */
					     sizeof(bad_randu_state_t), 
					     &bad_randu_set, 
					     &bad_randu_get } ;

const gsl_rng_type * gsl_rng_bad_randu = &bad_randu_type ;
