/* $Id$ */
/* This file was generated automatically from the template file
 * x-gen.xc (either that or it _is_ the template file)
 */

#include "XXX.h"             /* defines gsl_ran_XXX_ prototypes */
#include "gsl_ran.h"         /* defines gsl_ran_ prototypes */

inline unsigned long gsl_ran_random_wstate(void *vState) {
    return gsl_ran_XXX_random_wstate(vState);
}
inline double gsl_ran_uniform_wstate(void *vState) {
    return gsl_ran_XXX_uniform_wstate(vState);
}
inline double gsl_ran_max() {
    return gsl_ran_XXX_max();
}
inline void gsl_ran_seed_wstate(void *vState, int seed) {
    gsl_ran_XXX_seed_wstate(vState,seed);
}
inline unsigned long gsl_ran_random() {
    return gsl_ran_XXX_random();
}
inline double gsl_ran_uniform() {
    return gsl_ran_XXX_uniform();
}
inline void gsl_ran_seed(int seed) {
    gsl_ran_XXX_seed(seed);
}
inline void *gsl_ran_getRandomState() {
    return gsl_ran_XXX_getRandomState();
}
inline void gsl_ran_setRandomState(void *vState) {
    gsl_ran_XXX_setRandomState(vState);
}
inline void gsl_ran_printState(void *vState) {
     gsl_ran_XXX_printState(vState);
}




