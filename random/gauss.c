#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_ran.h"

inline double gsl_ran_gaussian_wstate(gsl_ran_gaussianRandomState *gState)
{
    /* if the state includes a value from the last call,
       then just return that */
    double x,y,rr;
    
    if (gState->ng == 1) {
        gState->ng = 0;
        return gState->g;
    }

    if (gState->randomState == (void *)0) {
        gState->randomState = gsl_ran_getRandomState();
    }

    do {
        /* choose x,y in uniform square */
        x = -1.0+2.0*gsl_ran_uniform_wstate(gState->randomState);
        y = -1.0+2.0*gsl_ran_uniform_wstate(gState->randomState);
        rr = x*x+y*y;
    } while( rr > 1.0 && !(x==0 && y==0) );

    rr = sqrt(-2.0*log(rr)/rr); /* Box-Muller transform */
    gState->g = x*rr;           /* save one of the random variates */
    gState->ng = 1;             /* indicate we have an extra one now */
    return y*rr;                /* return the other */
}
void gsl_ran_copyGaussState(gsl_ran_gaussianRandomState *t,
                            gsl_ran_gaussianRandomState *f)
{
    /* the copied state is not an element-by-element copy, because
     * we want to make a new copy of the randomState */
    /* Note also, we should probably put in some checking whether
     * these are valid pointers! (or at least non-null) */
    /* It might be more efficient to use gsl_ran_copyRandomState(),
     * but that is a 'private' function;  instead we save the static
     * state, write the 'from' state to the static state (setState),
     * copy the static state to the 'to' state (requires a malloc),
     * and then reset the static state to the saved tmp state.  Finally,
     * it is save to free the tmp state.
     */
    
    void *tmpRandomState;
    tmpRandomState = gsl_ran_getRandomState(); /* keep copy of original */
    gsl_ran_setRandomState(f->randomState);
    /* if its already been allocated, free it! */
    if (t->randomState != NULL)
        cfree((char *)t->randomState);
    t->randomState = gsl_ran_getRandomState();
    gsl_ran_setRandomState(tmpRandomState);    /* reset original */
    cfree((char *)tmpRandomState);
    
    t->ng = f->ng;
    t->g = f->g;
}

static gsl_ran_gaussianRandomState gstate = {0, 0.0, (void *)0};

double gsl_ran_gaussian() 
{
    return gsl_ran_gaussian_wstate(&gstate);
}
gsl_ran_gaussianRandomState *gsl_ran_getGaussState(void)
{
    gsl_ran_gaussianRandomState *theState;
    theState = (gsl_ran_gaussianRandomState *)calloc(1,sizeof(gsl_ran_gaussianRandomState));
    gsl_ran_copyGaussState(theState,&gstate);
    return theState;
}
void gsl_ran_setGaussState(gsl_ran_gaussianRandomState *theState)
{
    gsl_ran_copyGaussState(&gstate,theState);
}

    
