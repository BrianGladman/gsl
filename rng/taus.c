#include <config.h>
#include <stdlib.h>
#include <gsl_rng.h>

/* From:
   P. L'Ecuyer, "Maximally equidistributed combined Tausworthe
   generators" [***insert ref ***]
   This is available on the net from L'Ecuyer's home page.
   [*** insert URL ***]  */

unsigned long int taus_get (void * vstate);
void taus_set (void * state, unsigned int s);

typedef struct {
    unsigned long int s1, s2, s3;
} taus_state_t ;

unsigned long taus_get(void * vstate)
{
    taus_state_t * state = (taus_state_t *) vstate;

    /* GNU's gcc likes the 'UL' suffix to indicate unsigned long */

#define TAUSWORTHE(s,a,b,c,d) ((( s & c ) << d ) ^ ((( s << a) ^ s ) >> b))

    state->s1 = TAUSWORTHE(state->s1, 13, 19, 4294967294UL, 12);
    state->s2 = TAUSWORTHE(state->s2,  2, 25, 4294967288UL,  4);
    state->s3 = TAUSWORTHE(state->s3,  3, 11, 4294967280UL, 17);

    return (state->s1 ^ state->s2 ^ state->s3) ;
}

/* LCG is a "quick and dirty" (Press et al, p284) generator */ 
#define LCG(n) ((n)*8121+28411)%134456

void taus_set(void * vstate, unsigned int s)
{
  /* An entirely adhoc way of seeding!  L'Ecuyer suggests: s1,s2,s3 >=
     2,8,16, and says "In practice, it is better to take larger
     (random) initial seeds" */
  
  taus_state_t * state = (taus_state_t *) vstate;
    
  if (s == 0) s = 1;

  state->s1 = LCG(s);
  state->s2 = LCG(state->s1);
  state->s3 = LCG(state->s2);

  /* "warm it up" */
  taus_get (state);
  taus_get (state);
  taus_get (state);
  taus_get (state);
  taus_get (state);
  taus_get (state);
  return;
}

static const gsl_rng_type taus_type = { "gsl-taus",  /* name */
					4294967295UL,  /* RAND_MAX */
					sizeof(taus_state_t), 
					&taus_set, 
					&taus_get } ;

const gsl_rng_type * gsl_rng_taus = &taus_type ;


