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
void gsl_ran_copyGaussState(gsl_ran_gaussianRandomState *t, gsl_ran_gaussianRandomState *f)
{
    t->randomState = f->randomState;
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
    theState->randomState = gsl_ran_getRandomState();
    gsl_ran_copyState(theState->randomState,gstate.randomState);
    return theState;
}
void gsl_ran_setGaussState(gsl_ran_gaussianRandomState *theState)
{
    gsl_ran_copyGaussState(&gstate,theState);
    gstate.randomState = gsl_ran_getRandomState();
    gsl_ran_copyState(gstate.randomState,theState->randomState);
}

    
