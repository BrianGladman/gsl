#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_random.h"

inline double GSL_gauss_wstate(GSL_gaussRandomState *gState)
{
    /* if the state includes a value from the last call,
       then just return that */
    double x,y,rr;
    
    if (gState->ng == 1) {
        gState->ng = 0;
        return gState->g;
    }

    if (gState->randomState == (void *)0) {
        gState->randomState = GSL_getRandomState();
    }

    do {
        /* choose x,y in uniform square */
        x = -1.0+2.0*GSL_uniform_wstate(gState->randomState);
        y = -1.0+2.0*GSL_uniform_wstate(gState->randomState);
        rr = x*x+y*y;
    } while( rr > 1.0 && !(x==0 && y==0) );

    rr = sqrt(-2.0*log(rr)/rr); /* Box-Muller transform */
    gState->g = x*rr;           /* save one of the random variates */
    gState->ng = 1;             /* indicate we have an extra one now */
    return y*rr;                /* return the other */
}
void GSL_copyGaussState(GSL_gaussRandomState *t, GSL_gaussRandomState *f)
{
    t->randomState = f->randomState;
    t->ng = f->ng;
    t->g = f->g;
}

static GSL_gaussRandomState gstate = {0, 0.0, (void *)0};

double GSL_gauss() 
{
    return GSL_gauss_wstate(&gstate);
}
GSL_gaussRandomState *GSL_getGaussState(void)
{
    GSL_gaussRandomState *theState;
    theState = (GSL_gaussRandomState *)calloc(1,sizeof(GSL_gaussRandomState));
    GSL_copyGaussState(theState,&gstate);
    theState->randomState = GSL_getRandomState();
    GSL_copyRandomState(theState->randomState,gstate.randomState);
    return theState;
}
void GSL_setGaussState(GSL_gaussRandomState *theState)
{
    GSL_copyGaussState(&gstate,theState);
    gstate.randomState = GSL_getRandomState();
    GSL_copyRandomState(gstate.randomState,theState->randomState);
}

    
