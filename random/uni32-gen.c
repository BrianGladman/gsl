/* $Id$ */
/* File xxx-gen.xc is a template file
 *
 * This is either the template file xxx-gen.xc or else is generated
 * automatically from it.  If you want to make changes to this file,
 * make sure that the changes are to the template file xxx-gen.xc;
 * you can then rebuild the files uni32-switch.c, uni32-gen.c, and uni32.h by
 * running the script 'makesrc uni32' If you run the script with no
 * arguments, it makes all built-sources
 */


#include "uni32.h"             /* defines gsl_ran_uni32_ prototypes */
#include "gsl_ran.h"         /* defines gsl_ran_ prototypes */

inline unsigned long gsl_ran_random_wstate(void *vState) {
    return gsl_ran_uni32_random_wstate(vState);
}
inline double gsl_ran_uniform_wstate(void *vState) {
    return gsl_ran_uni32_uniform_wstate(vState);
}
inline double gsl_ran_max() {
    return gsl_ran_uni32_max();
}
inline void gsl_ran_seed_wstate(void *vState, int seed) {
    gsl_ran_uni32_seed_wstate(vState,seed);
}
inline unsigned long gsl_ran_random() {
    return gsl_ran_uni32_random();
}
inline double gsl_ran_uniform() {
    return gsl_ran_uni32_uniform();
}
inline void gsl_ran_seed(int seed) {
    gsl_ran_uni32_seed(seed);
}
inline void *gsl_ran_getRandomState() {
    return gsl_ran_uni32_getRandomState();
}
inline void gsl_ran_setRandomState(void *vState) {
    gsl_ran_uni32_setRandomState(vState);
}
inline void gsl_ran_printState(void *vState) {
     gsl_ran_uni32_printState(vState);
}




