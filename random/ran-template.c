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
#include "gsl_random.h"  /* also required: function prototypes, mostly */

/* Always define the GSL_randomState type */
typedef struct {
    /* Include whatever it takes to define the state */
    unsigned long s1,s2,s3;
} GSL_randomState;



/* inline usually saves some function call overhead */
inline unsigned long GSL_random_wstate(void *vState) {

    unsigned long r;
    
    GSL_randomState *theState;
    theState = (GSL_randomState *)vState;

    /* include code here: state variables are referred to as
     * theState->x1, etc... */

    return r;
}

inline double GSL_randmax(void)
{
    return XXXXXXXX;            /* usually some large constant */
}

void GSL_seed_wstate(int seed, void *vState)
{
    GSL_randomState *theState;
    theState = (GSL_randomState *)vState;

    /* Set the GSL_randomState to different random states according to
     * the different choices of seed.  (Even slightly) different seeds
     * should provide (extremely) different states.
     */

}

inline double GSL_uniform_wstate(void *theState)
{
    /* Some PRNG's work directly with double's, but most manipulate
     * integers, and then the uniform() routine looks something like
     * the following, where YYYYYYY is 1/GSL_randmax().
     */

    return (double)GSL_random(theState)*YYYYYYY;

}

void GSL_copyRandomState(void *tState, void *fState)
{
    GSL_randomState *toState, *fromState;
    toState   = (GSL_randomState *)tState;
    fromState = (GSL_randomState *)fState;

    /* This should be an element-by-element copy from the fromState
     * to the toState; eg */

    toState->s1 = fromState->s1; /* etc... */
}

/* Everything up to here has explicitly referred to the randomState in
 * the argument list.  What follows will make reference to the static
 * state that is kept here for convenience */

static GSL_randomState state = { /* include initialization */ };

/* if static initialization for the full GSL_randomState is not feasable, then
 * at least have an 'int inityet' component of the state which is initialized
 * to zero and then set to one by the GSL_seed_wstate() routine */

/* ranstate.c is boilerplate that defines:
   GSL_random(), GSL_uniform(), GSL_seed()
   void *GSL_getRandomState(),
   void GSL_setRandomState(void *)
   */

#include "ranstate.c"

/* ...that's all */
