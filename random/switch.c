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
 * expand 'taus' into the favorite 'xxx'
 */

void gsl_ran_use_default() {
    gsl_ran_use_taus();
}

/*
 * The following function is a template which the shell script will
 * expand into a set of funtions, one for each algorithm: mrg, cmrg, etc
 */
#include "taus.h"
void gsl_ran_use_taus() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("taus");
    A->random_wstate = gsl_ran_taus_random_wstate;
    A->uniform_wstate = gsl_ran_taus_uniform_wstate;
    A->max = gsl_ran_taus_max;
    A->random = gsl_ran_taus_random;
    A->uniform = gsl_ran_taus_uniform;
    A->seed = gsl_ran_taus_seed;
    A->getRandomState = gsl_ran_taus_getRandomState;
    A->setRandomState = gsl_ran_taus_setRandomState;
}

#include "mrg.h"
void gsl_ran_use_mrg() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("mrg");
    A->random_wstate = gsl_ran_mrg_random_wstate;
    A->uniform_wstate = gsl_ran_mrg_uniform_wstate;
    A->max = gsl_ran_mrg_max;
    A->random = gsl_ran_mrg_random;
    A->uniform = gsl_ran_mrg_uniform;
    A->seed = gsl_ran_mrg_seed;
    A->getRandomState = gsl_ran_mrg_getRandomState;
    A->setRandomState = gsl_ran_mrg_setRandomState;
}

#include "cmrg.h"
void gsl_ran_use_cmrg() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("cmrg");
    A->random_wstate = gsl_ran_cmrg_random_wstate;
    A->uniform_wstate = gsl_ran_cmrg_uniform_wstate;
    A->max = gsl_ran_cmrg_max;
    A->random = gsl_ran_cmrg_random;
    A->uniform = gsl_ran_cmrg_uniform;
    A->seed = gsl_ran_cmrg_seed;
    A->getRandomState = gsl_ran_cmrg_getRandomState;
    A->setRandomState = gsl_ran_cmrg_setRandomState;
}

#include "uni.h"
void gsl_ran_use_uni() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("uni");
    A->random_wstate = gsl_ran_uni_random_wstate;
    A->uniform_wstate = gsl_ran_uni_uniform_wstate;
    A->max = gsl_ran_uni_max;
    A->random = gsl_ran_uni_random;
    A->uniform = gsl_ran_uni_uniform;
    A->seed = gsl_ran_uni_seed;
    A->getRandomState = gsl_ran_uni_getRandomState;
    A->setRandomState = gsl_ran_uni_setRandomState;
}

#include "uni32.h"
void gsl_ran_use_uni32() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("uni32");
    A->random_wstate = gsl_ran_uni32_random_wstate;
    A->uniform_wstate = gsl_ran_uni32_uniform_wstate;
    A->max = gsl_ran_uni32_max;
    A->random = gsl_ran_uni32_random;
    A->uniform = gsl_ran_uni32_uniform;
    A->seed = gsl_ran_uni32_seed;
    A->getRandomState = gsl_ran_uni32_getRandomState;
    A->setRandomState = gsl_ran_uni32_setRandomState;
}

#include "zuf.h"
void gsl_ran_use_zuf() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("zuf");
    A->random_wstate = gsl_ran_zuf_random_wstate;
    A->uniform_wstate = gsl_ran_zuf_uniform_wstate;
    A->max = gsl_ran_zuf_max;
    A->random = gsl_ran_zuf_random;
    A->uniform = gsl_ran_zuf_uniform;
    A->seed = gsl_ran_zuf_seed;
    A->getRandomState = gsl_ran_zuf_getRandomState;
    A->setRandomState = gsl_ran_zuf_setRandomState;
}

#include "rand.h"
void gsl_ran_use_rand() {
    if (A==NULL) gsl_ran_newAlgorithm();

    A->name = strdup("rand");
    A->random_wstate = gsl_ran_rand_random_wstate;
    A->uniform_wstate = gsl_ran_rand_uniform_wstate;
    A->max = gsl_ran_rand_max;
    A->random = gsl_ran_rand_random;
    A->uniform = gsl_ran_rand_uniform;
    A->seed = gsl_ran_rand_seed;
    A->getRandomState = gsl_ran_rand_getRandomState;
    A->setRandomState = gsl_ran_rand_setRandomState;
}

