#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* MINSTD is Park and Miller's minimal standard generator (i.e. it's
   not particularly good).

   The sequence is
   
         x_{n+1} = (a x_n) mod m

   with a = 16807 and m = 2^31 - 1 = 2147483647. The seed specifies
   the initial value, x_0.  

   From: Park and Miller, "Random Number Generators: Good ones are
   hard to find" Communications of the ACM, October 1988, Volume 31,
   No 10, pages 1192-1201. */

unsigned long int minstd_get (void * vstate);
void minstd_set (void * state, unsigned int s);
void minstd_set_with_state (void * vstate, const void * vinit_state,
			 unsigned int s);

static const int m = 2147483647, a = 16807, q = 127773, r = 2836 ;

typedef struct {
  unsigned int x ;
} minstd_state_t ;

unsigned long int minstd_get (void *vstate)
{
    minstd_state_t * state = (minstd_state_t *)vstate;

    const unsigned int x = state->x ;
    
    const int h  = x / q;    
    const int t = a * (x - h * q) - h * r;

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


static const minstd_state_t init_state = {
    1
};

void 
minstd_set(void * vstate, unsigned int s)
{
  minstd_state_t * state = (minstd_state_t *) vstate;
  
  if (s == 0) s = 1;

  state->x = s ;
  
  return ;
}

static const gsl_rng_type minstd_type = { "minstd",  /* name */
					  2147483647,  /* RAND_MAX */
					  sizeof(minstd_state_t), 
					  &minstd_set, 
					  &minstd_get } ;

const gsl_rng_type * gsl_rng_minstd = &minstd_type ; 
