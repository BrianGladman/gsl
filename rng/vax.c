#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* This is the old vax generator MTH$RANDOM. The sequence is,
   
         x_{n+1} = (a x_n + c) mod m

   with a = 69069, c = 1 and m = 2^32. The seed specifies the initial
   value, x_1.

   The theoretical value of x_{10001} is 3051034865.

   The period of this generator is 2^32. */

unsigned long int vax_get (void * vstate);
void vax_set (void * state, unsigned long int s);
void vax_set_with_state (void * vstate, const void * vinit_state,
			 unsigned long int s);

static const long int a = 69069 ;
static const long int c = 1 ;
static const unsigned long int m = 4294967295UL ;
static const long int q = 62183 ;
static const long int r = 49669 ;

typedef struct {
  unsigned long int x;
} vax_state_t ;

unsigned long int vax_get (void *vstate)
{
    vax_state_t * state = (vax_state_t *)vstate;

    const unsigned long int x = state->x ;

    const long int h = x / q ;
    const long int t = a * (x - h * q) - h * r + c;

    if (t < 0) 
      {
	state->x = t + m + 1;
      }
    else
      {
	state->x = t ;
      }
    
    return state->x;
}

void vax_set(void * vstate, unsigned long int s)
{
  vax_state_t * state = (vax_state_t *) vstate;
  
  state->x = s ;

  return;
}

static const gsl_rng_type vax_type = { "vax",  /* name */
				       4294967295UL,  /* RAND_MAX */
				       sizeof(vax_state_t), 
				       &vax_set, 
				       &vax_get } ;

const gsl_rng_type * gsl_rng_vax = &vax_type ;
