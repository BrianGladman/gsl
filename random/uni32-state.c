/* File xxx-state.xc is a template file
 *
 * This is either the template file xxx-state.xc or else is generated
 * automatically from it.  If you want to make changes to this file,
 * make sure that the changes are to the template file xxx-state.xc;
 * you can then rebuild the files XXX-state.c, XXX-gen.c, and XXX.h by
 * running the script 'makesrc XXX'.  If you run the script with as
 * 'makesrc --all' then it makes all built-sources
 */
/*
 * This file defines functions gsl_ran_uni32_
 *   seed(),random(),uniform(),copyState(),getState(),setState()
 * in terms of the routines defined by the human programmer in XXX.c
 *
 * It is #include'd in file XXX.c
 */

void gsl_ran_uni32_seed(int s)
{
    gsl_ran_uni32_seed_wstate((void *)&state,s);
}
inline void gsl_ran_uni32_printState(void *s)
{
    gsl_ran_uni32_printState_p((gsl_ran_uni32_randomState *)s);
}
inline unsigned long gsl_ran_uni32_random(void)
{
    return gsl_ran_uni32_random_wstate((void *)&state);
}
static const double gsl_ran_uni32_invRANDMAX = 1.0/(1.0+gsl_ran_uni32_RANDMAX);
double gsl_ran_uni32_max() { return (double)gsl_ran_uni32_RANDMAX; }
inline double gsl_ran_uni32_uniform_wstate(void *s)
{
    /* Returned value is guaranteed to be strictly between zero and one */
    double Z;
    Z = gsl_ran_uni32_random_wstate(s);
    if (Z==0) Z=gsl_ran_uni32_RANDMAX;
    return Z*gsl_ran_uni32_invRANDMAX;
}
inline double gsl_ran_uni32_uniform(void)
{
    return gsl_ran_uni32_uniform_wstate((void *)&state);
}

/* get/set randomState */
/* get/set randomState */
static inline void gsl_ran_uni32_copyState(void *tState,void *fState)
{
    bcopy((char *)fState,(char*)tState,sizeof(gsl_ran_uni32_randomState));
}

void *gsl_ran_uni32_getRandomState(void)
{
    gsl_ran_uni32_randomState *theState;
    theState = (gsl_ran_uni32_randomState *)calloc(1,sizeof(gsl_ran_uni32_randomState));
    gsl_ran_uni32_copyState(theState,(void *)&state);
    return (void *)theState;
}
void gsl_ran_uni32_setRandomState(void *vState)
{
    gsl_ran_uni32_copyState((void *)&state,(gsl_ran_uni32_randomState *)vState);
}


