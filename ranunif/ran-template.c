/* $Id$ */
/* This is not real code, but a template for adding new PRNG's (pseudo
 * random number generators) consistent with the GSL interface.
 *
 * In this space, you should copiously document your new PRNG.
 * In particular, where did you get it.  Where was it published,
 * etc.
 */

/* Standard include's: */
#include <stdio.h>              /* you may not need these two. */
#include <math.h>

#include <stdlib.h>      /* you'll need this one for calloc() which
                          * is how new randomState's are allocated */
#include "gsl_ran.h"  /* also required: function prototypes, mostly */

/* Always define the gsl_ran_randomState type */
typedef struct {
    /* Include whatever it takes to define the state */
    unsigned long s1,s2,s3;
} gsl_ran_randomState;



/* inline usually saves some function call overhead */
inline unsigned long gsl_ran_random_wstate(void *vState) {

    unsigned long r;
    
    gsl_ran_randomState *theState;
    theState = (gsl_ran_randomState *)vState;

    /* include code here: state variables are referred to as
     * theState->x1, etc... */

    return r;
}

inline double gsl_ran_max(void)
{
    return XXXXXXXX;            /* usually some large constant */
}

void gsl_ran_seed_wstate(int seed, void *vState)
{
    gsl_ran_randomState *theState;
    theState = (gsl_ran_randomState *)vState;

    /* Set the gsl_ran_randomState to different random states according to
     * the different choices of seed.  (Even slightly) different seeds
     * should provide (extremely) different states.
     */

}

inline double gsl_ran_uniform_wstate(void *theState)
{
    /* Some PRNG's work directly with double's, but most manipulate
     * integers, and then the uniform() routine looks something like
     * the following, where YYYYYYY is 1/gsl_ran_max().
     */

    return (double)gsl_ran_random(theState)*YYYYYYY;

}

void gsl_ran_copyState(void *tState, void *fState)
{
    gsl_ran_randomState *toState, *fromState;
    toState   = (gsl_ran_randomState *)tState;
    fromState = (gsl_ran_randomState *)fState;

    /* This should be an element-by-element copy from the fromState
     * to the toState; eg */

    toState->s1 = fromState->s1; /* etc... */
}

/* Everything up to here has explicitly referred to the randomState in
 * the argument list.  What follows will make reference to the static
 * state that is kept here for convenience */

static gsl_ran_randomState state = { /* include initialization */ };

/* if static initialization for the full gsl_ran_randomState is not feasable, then
 * at least have an 'int inityet' component of the state which is initialized
 * to zero and then set to one by the gsl_ran_seed_wstate() routine */

/* ranstate.c is boilerplate that defines:
   gsl_ran_random(), gsl_ran_uniform(), gsl_ran_seed()
   void *gsl_ran_getRandomState(),
   void gsl_ran_setRandomState(void *)
   */

#include "ranstate.c"

/* ...that's all */
