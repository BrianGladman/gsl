/* $Id$ */
/**
  From:
  P. L'Ecuyer, "Combined Multiple Recursive Random Number Generators,"
  to appear in Operations Research, 1996.
  (Preprint obtained as file compmrg.ps from L'Ecuyer's web page.)
  **/

#include <stdlib.h>
#include "gsl_random.h"

const int m1 = 2147483647, m2 = 2145483479;
int a12=63308,   q12=33921, r12=12979,
    a13=-183326, q13=11714, r13=2883,
    a21=86098,   q21=24919, r21=7417,
    a23=-539608, q23=3976,  r23=2071;


typedef struct {
    long x10,x11,x12;           /* first component */
    long x20,x21,x22;           /* second component */
} GSL_randomState;

const double Invmp1 = 4.656612873077393e-10;

#define POSITIVE(x,m) if (x<0) x += m

inline unsigned long GSL_random_wstate(void *vState)
{
    int h,p12,p13,p21,p23;
    GSL_randomState *theState;
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
inline double GSL_uniform_wstate(void *vState)
{
    int Z;
    Z = GSL_random_wstate(vState);
    if (Z==0) Z=m1;
    return(Z*Invmp1);
}

double GSL_randmax() { return (double)m1; }

void GSL_seed_wstate(void *vState, int s)
{
    GSL_randomState *theState;
    theState = (GSL_randomState *)vState;

    /* An entirely adhoc way of seeding! This does not come
       from L'Ecuyer et al */
    theState->x10 = s;
    theState->x11 = 128+theState->x10;
    theState->x12 = 128+theState->x11;
    theState->x20 = 128+theState->x12;
    theState->x21 = 128+theState->x20;
    theState->x22 = 128+theState->x21;
    return;
}

void GSL_copyRandomState(void *tState,
                         void *fState)
{
    GSL_randomState *toState, *fromState;
    toState   = (GSL_randomState *)tState;
    fromState = (GSL_randomState *)fState;

    toState->x10 = fromState->x10;
    toState->x11 = fromState->x11;
    toState->x12 = fromState->x12;
    toState->x20 = fromState->x20;
    toState->x21 = fromState->x21;
    toState->x22 = fromState->x22;
}
    
/* -- From here on in, it's boilerplate! -- */

static GSL_randomState state;
#include "ranstate.c"

    
    
    
