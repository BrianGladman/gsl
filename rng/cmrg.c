/**
  From:
  P. L'Ecuyer, "Combined Multiple Recursive Random Number Generators,"
  to appear in Operations Research, 1996.
  (Preprint obtained as file compmrg.ps from L'Ecuyer's web page.)
  **/

#include <stdlib.h>
#include "gsl_ran.h"

static const int m1 = 2147483647, m2 = 2145483479;
static int a12=63308,   q12=33921, r12=12979,
           a13=-183326, q13=11714, r13=2883,
           a21=86098,   q21=24919, r21=7417,
           a23=-539608, q23=3976,  r23=2071;

#define gsl_ran_cmrg_RANDMAX 2147483647 /* m1 */

typedef struct {
    long x10,x11,x12;           /* first component */
    long x20,x21,x22;           /* second component */
} gsl_ran_cmrg_randomState;

static void
gsl_ran_cmrg_printState_p(gsl_ran_cmrg_randomState *s)
{
    printf("%ldL, %ldL, %ldL,\n%ldL, %ldL, %ldL\n",
	   s->x10,s->x11,s->x12,
	   s->x20,s->x21,s->x22);
}


#define POSITIVE(x,m) if (x<0) x += m

inline unsigned long gsl_ran_cmrg_random_wstate(void *vState)
{
    int h,p12,p13,p21,p23;
    gsl_ran_cmrg_randomState *theState;
    theState = (void *)vState;

    /* Component 1 */
    h = theState->x10 / q13; p13 = -a13  *(theState->x10-h*q13) - h*r13;
    h = theState->x11 / q12; p12 =  a12  *(theState->x11-h*q12) - h*r12;
    POSITIVE(p13,m1);
    POSITIVE(p12,m1);
    theState->x10 = theState->x11;
    theState->x11 = theState->x12;
    theState->x12 = p12-p13;
    POSITIVE(theState->x12,m1);
    
    /* Component 2 */
    h = theState->x20 / q23; p23 = -a23 * (theState->x20-h*q23)  - h*r23;
    h = theState->x22 / q21; p21 =  a21 * (theState->x22-h*q21)  - h*r21;
    POSITIVE(p23,m2);
    POSITIVE(p21,m2);
    theState->x20 = theState->x21;
    theState->x21 = theState->x22;
    theState->x22 = p21-p23;
    POSITIVE(theState->x22,m2);
    /* Combination */
    if(theState->x12 < theState->x22)
        return (theState->x12-theState->x22+m1);
    else
        return (theState->x12-theState->x22);
}

#define LCG(n) ((n)*8121+28411)%134456
void gsl_ran_cmrg_seed_wstate(void *vState, int s)
{
    /* An entirely adhoc way of seeding! This does **not** come
       from L'Ecuyer et al */
    gsl_ran_cmrg_randomState *rState;
    rState = (gsl_ran_cmrg_randomState *)vState;

    s = (s<0 ? -s : s);
    if (s==0) s=1;
    rState->x10 = LCG(s);
    rState->x11 = LCG(rState->x10);
    rState->x12 = LCG(rState->x11);
    rState->x20 = LCG(rState->x12);
    rState->x21 = LCG(rState->x20);
    rState->x22 = LCG(rState->x21);

    /* "warm it up" */
    gsl_ran_cmrg_random_wstate(rState);
    gsl_ran_cmrg_random_wstate(rState);
    gsl_ran_cmrg_random_wstate(rState);
    gsl_ran_cmrg_random_wstate(rState);
    gsl_ran_cmrg_random_wstate(rState);
    gsl_ran_cmrg_random_wstate(rState);
    gsl_ran_cmrg_random_wstate(rState);
    return;
}

static gsl_ran_cmrg_randomState state = {
    511515612L, 1645048169L, 1860274777L,
     55882945L, 1225790668L, 2055528708L
};
#include "cmrg-state.c"

    
    
    
