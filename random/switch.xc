/* $Id$ */
/* Switch between various algorithms at run-time, using the function
   gsl_ran_use_xxx() to start using xxx.  The function gsl_ran_name()
   returns the name of the algorithm currently in use.
   */

#include "gsl_ran.h"         /* defines gsl_ran_ prototypes */

#define MAXALGNAMELEN 10        /* only short algorithm names allowed */
typedef struct {
    char name[MAXALGNAMELEN];
    unsigned long (*random_wstate)(void *);
    double (*uniform_wstate)(void *);
    double (*max)();
    void (*seed_wstate)(void *,int);
    unsigned long (*random)();
    double (*uniform)();
    void (*seed)();
    void (*copyState)(void *,void *);
    void *(*getRandomState)();
    void (*setRandomState)(void *);    
} AlgorithmSwitch;

AlgorithmSwitch A;

char *gsl_ran_name() {
    return A.name;
}

inline unsigned long gsl_ran_random_wstate(void *vState) {
    return A.random_wstate(vState);
}
inline double gsl_ran_uniform_wstate(void *vState) {
    return A.uniform_wstate(vState);
}
inline double gsl_ran_max() {
    return A.max();
}
inline void gsl_ran_seed_wstate(void *vState, int seed) {
    A.seed_wstate(vState,seed);
}
inline unsigned long gsl_ran_random() {
    return A.random();
}
inline double gsl_ran_uniform() {
    return A.uniform();
}
inline void gsl_ran_seed(int seed) {
    A.seed(seed);
}
inline void gsl_ran_copyState(void *t, void *f) {
    A.copyState(t,f);
}
inline void *gsl_ran_getRandomState() {
    return A.getRandomState();
}
inline void gsl_ran_setRandomState(void *vState) {
    A.setRandomState(vState);
}

/*
 * The following function is a template which the shell script will
 * expand into a set of funtions, one for each algorithm: mrg, cmrg, etc
 */
#include "xxx.h"
void gsl_ran_use_xxx() {
    strncpy(A.name,"xxx",MAXALGNAMELEN);
    A.random_wstate = gsl_ran_xxx_random_wstate;
    A.uniform_wstate = gsl_ran_xxx_uniform_wstate;
    A.max = gsl_ran_xxx_max;
    A.random = gsl_ran_xxx_random;
    A.uniform = gsl_ran_xxx_uniform;
    A.seed = gsl_ran_xxx_seed;
    A.copyState = gsl_ran_xxx_copyState;
    A.getRandomState = gsl_ran_xxx_getRandomState;
    A.setRandomState = gsl_ran_xxx_setRandomState;
}

