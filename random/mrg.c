/* $Id$ */
/**
  From:
  P. L'Ecuyer, F. Blouin, and R. Coutre, "A search for good multiple
  recursive random number generators, ACM Transactions on Modeling
  and Computer Simulation 3, 87-98 (1998).
  **/

#include <stdlib.h>
#include "mrg.h"

static const long m = 2147483647,
    a1 = 107374182,   q1 = 20,      r1 = 7,
    a5 = 104480,      q5 = 20554,   r5 = 1727;

static const double  Invmp1 = 4.656612873077393e-10; /* = 1/m */

typedef struct {
    long x1,x2,x3,x4,x5;
} gsl_ran_mrg_randomState;

static gsl_ran_mrg_randomState state = { 12345, 23456, 34567, 45678, 56789 };

inline unsigned long gsl_ran_mrg_random_wstate(void *vState)
{
    long h, p1, p5;
    gsl_ran_mrg_randomState *theState;
    theState = (void *)vState;
    
    h  = theState->x5 / q5;    p5 = a5 * (theState->x5 - h * q5) - h * r5;
    theState->x5 = theState->x4;
    theState->x4 = theState->x3;
    theState->x3 = theState->x2;
    theState->x2 = theState->x1;
    h  = theState->x1 / q1;
    p1 = a1 * (theState->x1 - h * q1) - h * r1;
    if (p1 < 0) p1 += m;
    if (p5 > 0) p5 -= m;
    theState->x1 = p1 + p5;
    if (theState->x1 < 0) theState->x1 += m;
    return theState->x1;
}

inline double gsl_ran_mrg_uniform_wstate(void *vState)
{
    long Z;
    Z = gsl_ran_mrg_random_wstate(vState);
    if (Z == 0) Z = m;
    return (Z * Invmp1);
}

inline double gsl_ran_mrg_max(void)
{
    return (double)m;
}
void gsl_ran_mrg_seed_wstate(void *vState, int s)
{
    
    /* An entirely adhoc way of seeding! This does not come
       from L'Ecuyer et al */
    s -= 1;
    ((gsl_ran_mrg_randomState *)vState)->x1 = 12345 + s;
    ((gsl_ran_mrg_randomState *)vState)->x2 = 23456 + 2*s;
    ((gsl_ran_mrg_randomState *)vState)->x3 = 34567 + 3*s;
    ((gsl_ran_mrg_randomState *)vState)->x4 = 45678 + 5*s;
    ((gsl_ran_mrg_randomState *)vState)->x5 = 56789 + 7*s;
    return;
}


/* get/set randomState */
void gsl_ran_mrg_copyState(void *tState,
                         void *fState)
{
    gsl_ran_mrg_randomState *toState, *fromState;
    toState   = (gsl_ran_mrg_randomState *)tState;
    fromState = (gsl_ran_mrg_randomState *)fState;

    toState->x1 = fromState->x1;
    toState->x2 = fromState->x2;
    toState->x3 = fromState->x3;
    toState->x4 = fromState->x4;
    toState->x5 = fromState->x5;

}

static gsl_ran_mrg_randomState state;
#include "mrg-state.c"

    
    
