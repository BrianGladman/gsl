#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* From:
   P. L'Ecuyer, F. Blouin, and R. Coutre, "A search for good multiple
   recursive random number generators, ACM Transactions on Modeling
   and Computer Simulation 3, 87-98 (1998).  */

unsigned long int mrg_get (void * vstate);
void mrg_set (void * state, unsigned int s);
void mrg_set_with_state (void * vstate, const void * vinit_state,
			 unsigned int s);

static const int m = 2147483647,
  a1 = 107374182,   q1 = 20,      r1 = 7,
  a5 = 104480,      q5 = 20554,   r5 = 1727;


typedef struct {
  long int x1, x2, x3, x4, x5;
} mrg_state_t ;

unsigned long int mrg_get (void *vstate)
{
    long int  h, p1, p5;

    mrg_state_t * state = (mrg_state_t *)vstate;
    
    h  = state->x5 / q5;    
    p5 = a5 * (state->x5 - h * q5) - h * r5;
    state->x5 = state->x4;
    state->x4 = state->x3;
    state->x3 = state->x2;
    state->x2 = state->x1;
    h  = state->x1 / q1;
    p1 = a1 * (state->x1 - h * q1) - h * r1;
    if (p1 < 0) p1 += m;
    if (p5 > 0) p5 -= m;
    state->x1 = p1 + p5;
    if (state->x1 < 0) state->x1 += m;

    return state->x1;
}


static const mrg_state_t init_state = {
    1436981205L, 651435938L, 1374493895L, 1070522304L, 1168302460L
};

void mrg_set(void * state, unsigned int s)
{
  mrg_set_with_state (state, &init_state, s) ;
}

#define LCG(n) ((n)*8121+28411)%134456
void mrg_set_with_state (void * vstate, const void * vinit_state, unsigned int s)
{
  /* An entirely adhoc way of seeding! This does **not** come
     from L'Ecuyer et al */
  
  mrg_state_t * state = (mrg_state_t *) vstate;
  
  *state = *(const mrg_state_t *) vinit_state ;
  
  if (s == 0) s = 1;
  
  state->x1 = LCG(s);
  state->x2 = LCG(state->x1);
  state->x3 = LCG(state->x2);
  state->x4 = LCG(state->x3);
  state->x5 = LCG(state->x4);

  /* "warm it up" with at least 5 calls to go through
     all the x values */

  mrg_get(state);
  mrg_get(state);
  mrg_get(state);
  mrg_get(state);
  mrg_get(state);
  mrg_get(state);

  return;
}

static const gsl_rng_type mrg_type = { "gsl-mrg",  /* name */
					2147483647,  /* RAND_MAX */
					sizeof(mrg_state_t), 
					&mrg_set, 
					&mrg_get } ;

const gsl_rng_type * gsl_rng_mrg = &mrg_type ; 
