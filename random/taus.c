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

#define TAUSWORTHE(s,a,b,c,d) ((( s & c ) << d ) ^ ((( s << a) ^ s ) >> b))

typedef struct {
    unsigned long s1,s2,s3;
} gsl_ran_taus_randomState;

inline unsigned long gsl_ran_taus_random_wstate(void *vState)
{
    gsl_ran_taus_randomState *newState;
    newState = (gsl_ran_taus_randomState *)vState;
    /* GNU's gcc likes the 'UL' suffix to indicate unsigned long */
    newState->s1 = TAUSWORTHE(newState->s1, 13, 19, 4294967294UL, 12);
    newState->s2 = TAUSWORTHE(newState->s2,  2, 25, 4294967288UL,  4);
    newState->s3 = TAUSWORTHE(newState->s3,  3, 11, 4294967280UL, 17);

    return newState->s1 ^ newState->s2 ^ newState->s3;
}
inline double gsl_ran_taus_uniform_wstate(void *vState)
{
    return (double)gsl_ran_taus_random_wstate(vState)*2.3283064365e-10;
}
double gsl_ran_taus_max()
{
    return 4294967296.0;          /* = 2^32 */
}

void gsl_ran_taus_seed_wstate(void *vState,int s)
{
    /* An entirely adhoc way of seeding!
       L'Ecuyer suggests: s1,s2,s3 >= 2,8,16, and says
       "In practice, it is better to take larger (random) initial seeds"
       */
    s -= 1;
    ((gsl_ran_taus_randomState *)vState)->s1 = 12345 + s;
    ((gsl_ran_taus_randomState *)vState)->s2 = 23456 + 2*s;
    ((gsl_ran_taus_randomState *)vState)->s3 = 34567 + 3*s;
}



/* get/set randomState */
void gsl_ran_taus_copyState(void *tState,
                         void *fState)
{
    gsl_ran_taus_randomState *toState, *fromState;
    toState   = (gsl_ran_taus_randomState *)tState;
    fromState = (gsl_ran_taus_randomState *)fState;

    toState->s1 = fromState->s1;
    toState->s2 = fromState->s2;
    toState->s3 = fromState->s3;

}
    
    
static gsl_ran_taus_randomState state = { 12345, 54321, 98765 };
#include "taus-state.c"


