/* $Id$ */
/**
  From:
  P. L'Ecuyer, "Maximally equidistributed combined Tausworthe
  generators" [***insert ref ***]
  This is available on the net from L'Ecuyer's home page.
  [*** insert URL ***]
  **/

#include <stdlib.h>             /* calloc() */
#include "gsl_ran.h"
#define gsl_ran_taus_RANDMAX 4294967295.0 /* = 2^32  */

#define TAUSWORTHE(s,a,b,c,d) ((( s & c ) << d ) ^ ((( s << a) ^ s ) >> b))

typedef struct {
    unsigned long s1,s2,s3;
} gsl_ran_taus_randomState;

static void
gsl_ran_taus_printState_p(gsl_ran_taus_randomState *s)
{
    printf("%luUL, %luUL, %luUL\n",
	   s->s1,s->s2,s->s3);
}

unsigned long gsl_ran_taus_random_wstate(void *vState)
{
    gsl_ran_taus_randomState *newState;
    newState = (gsl_ran_taus_randomState *)vState;
    /* GNU's gcc likes the 'UL' suffix to indicate unsigned long */
    newState->s1 = TAUSWORTHE(newState->s1, 13, 19, 4294967294UL, 12);
    newState->s2 = TAUSWORTHE(newState->s2,  2, 25, 4294967288UL,  4);
    newState->s3 = TAUSWORTHE(newState->s3,  3, 11, 4294967280UL, 17);

    return newState->s1 ^ newState->s2 ^ newState->s3;
}

/* LCG is a "quick and dirty" (Press et al, p284) generator 
 */ 
#define LCG(n) ((n)*8121+28411)%134456
void gsl_ran_taus_seed_wstate(void *vState, int s)
{
    /* An entirely adhoc way of seeding!
       L'Ecuyer suggests: s1,s2,s3 >= 2,8,16, and says
       "In practice, it is better to take larger (random) initial seeds"
       */
    gsl_ran_taus_randomState *rState;
    rState = (gsl_ran_taus_randomState *)vState;

    s = (s<0 ? -s : s);
    if (s==0) s=1;
    rState->s1 = LCG(s);
    rState->s2 = LCG(rState->s1);
    rState->s3 = LCG(rState->s2);

    /* "warm it up" */
    gsl_ran_taus_random_wstate(rState);
    gsl_ran_taus_random_wstate(rState);
    gsl_ran_taus_random_wstate(rState);
    gsl_ran_taus_random_wstate(rState);
    gsl_ran_taus_random_wstate(rState);
    gsl_ran_taus_random_wstate(rState);
    return;
}

static gsl_ran_taus_randomState state = { 
    956008634UL, 2013275612UL, 4211549042UL
};
#include "taus-state.c"


