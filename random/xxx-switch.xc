/* $Id$ */
/* Switch between various algorithms at run-time, using the function
   gsl_ran_use_xxx() to start using xxx.  The function gsl_ran_name()
   returns the name of the algorithm currently in use.
   */

#include <stdio.h>              /* defines NULL */
#include <string.h>             /* for strdup() */
#include "gsl_ran.h"            /* defines gsl_ran_ prototypes */
#include "gsl_ran_switch.h"

typedef struct {
    char *name;
    unsigned long (*random_wstate)(void *);
    double (*uniform_wstate)(void *);
    double (*max)();
    void (*seed_wstate)(void *,int);
    unsigned long (*random)();
    double (*uniform)();
    void (*seed)();
    void *(*getRandomState)();
    void (*setRandomState)(void *);    
} AlgorithmSwitch;

AlgorithmSwitch *A=NULL;

inline void
gsl_ran_newAlgorithm()
{
    if (A == NULL) {
        A = (AlgorithmSwitch *)malloc(sizeof(AlgorithmSwitch));
    }
    gsl_ran_use_default();
}

char *gsl_ran_name() {
    if (A==NULL) gsl_ran_newAlgorithm();
    return A->name;
}

inline unsigned long gsl_ran_random_wstate(void *vState) {
    if (A==NULL) gsl_ran_newAlgorithm();
    return A->random_wstate(vState);
}
inline double gsl_ran_uniform_wstate(void *vState) {
    if (A==NULL) gsl_ran_newAlgorithm();
    return A->uniform_wstate(vState);
}
inline double gsl_ran_max() {
    if (A==NULL) gsl_ran_newAlgorithm();
    return A->max();
}
inline void gsl_ran_seed_wstate(void *vState, int seed) {
    if (A==NULL) gsl_ran_newAlgorithm();
    A->seed_wstate(vState,seed);
}
inline unsigned long gsl_ran_random() {
    if (A==NULL) gsl_ran_newAlgorithm();
    return A->random();
}
inline double gsl_ran_uniform() {
    if (A==NULL) gsl_ran_newAlgorithm();
    return A->uniform();
}
inline void gsl_ran_seed(int seed) {
    if (A==NULL) gsl_ran_newAlgorithm();
    A->seed(seed);
}
inline void *gsl_ran_getRandomState() {
    if (A==NULL) gsl_ran_newAlgorithm();
    return A->getRandomState();
}
inline void gsl_ran_setRandomState(void *vState) {
    if (A==NULL) gsl_ran_newAlgorithm();
    A->setRandomState(vState);
}

/*
 * This function is a template in which the shell script will
 * expand 'YYY' into the favorite 'xxx'
 */

void gsl_ran_use_default() {
    gsl_ran_use_YYY();
}

/*
 * The following function is a template which the shell script will
 * expand into a set of funtions, one for each algorithm: mrg, cmrg, etc
 */
#include "XXX.h"
void gsl_ran_use_XXX() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("XXX");
    A->random_wstate = gsl_ran_XXX_random_wstate;
    A->uniform_wstate = gsl_ran_XXX_uniform_wstate;
    A->max = gsl_ran_XXX_max;
    A->random = gsl_ran_XXX_random;
    A->uniform = gsl_ran_XXX_uniform;
    A->seed = gsl_ran_XXX_seed;
    A->getRandomState = gsl_ran_XXX_getRandomState;
    A->setRandomState = gsl_ran_XXX_setRandomState;
}

