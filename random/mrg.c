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

#define gsl_ran_mrg_RANDMAX 2147483647

typedef struct {
    long x1,x2,x3,x4,x5;
} gsl_ran_mrg_randomState;

static void
gsl_ran_mrg_printState_p(gsl_ran_mrg_randomState *s)
{
  printf("%ldL, %ldL, %ldL, %ldL, %ldL\n",
	   s->x1,s->x2,s->x3,s->x4,s->x5);
}

unsigned long gsl_ran_mrg_random_wstate(void *vState)
{
    long h, p1, p5;
    gsl_ran_mrg_randomState *theState;
    theState = (void *)vState;
    
    h  = theState->x5 / q5;    
    p5 = a5 * (theState->x5 - h * q5) - h * r5;
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

/* LCG is a "quick and dirty" (Press et al, p284) generator 
 */ 
#define LCG(n) ((n)*8121+28411)%134456
void gsl_ran_mrg_seed_wstate(void *vState, int s)
{
    /* An entirely adhoc way of seeding! This does **not** come
       from L'Ecuyer et al */
    gsl_ran_mrg_randomState *rState;
    rState = (gsl_ran_mrg_randomState *)vState;

    s = (s<0 ? -s : s);
    if (s==0) s=1;
    rState->x1 = LCG(s);
    rState->x2 = LCG(rState->x1);
    rState->x3 = LCG(rState->x2);
    rState->x4 = LCG(rState->x3);
    rState->x5 = LCG(rState->x4);
    /* "warm it up" with at least 5 calls to go through
     * all the x values */
    gsl_ran_mrg_random_wstate(rState);
    gsl_ran_mrg_random_wstate(rState);
    gsl_ran_mrg_random_wstate(rState);
    gsl_ran_mrg_random_wstate(rState);
    gsl_ran_mrg_random_wstate(rState);
    gsl_ran_mrg_random_wstate(rState);
    return;
}


static gsl_ran_mrg_randomState state = {
    1436981205L, 651435938L, 1374493895L, 1070522304L, 1168302460L
};



#include "mrg-state.c"

    
    
