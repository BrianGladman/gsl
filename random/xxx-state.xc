/* $Id$ */
/* File xxx-state.xc is a template file
 * This is either file xxx-state.xc or else is generated automatically 
 * from the template file xxx-state.xc
 *
 * It defines functions gsl_ran_XXX_
 *   seed(),random(),uniform(),copyState(),getState(),setState()
 * in terms of the routines defined by the human programmer in XXX.c
 *
 * It is #include'd in file XXX.c
 */

void gsl_ran_XXX_seed(int s)
{
    gsl_ran_XXX_seed_wstate((void *)&state,s);
}
inline void gsl_ran_XXX_printState(void *s)
{
    gsl_ran_XXX_printState_p((gsl_ran_XXX_randomState *)s);
}
inline unsigned long gsl_ran_XXX_random(void)
{
    return gsl_ran_XXX_random_wstate((void *)&state);
}
static const double gsl_ran_XXX_invRANDMAX = 1.0/(1.0+gsl_ran_XXX_RANDMAX);
double gsl_ran_XXX_max() { return (double)gsl_ran_XXX_RANDMAX; }
inline double gsl_ran_XXX_uniform_wstate(void *s)
{
    /* Returned value is guaranteed to be strictly between zero and one */
    double Z;
    Z = gsl_ran_XXX_random_wstate(s);
    if (Z==0) Z=gsl_ran_XXX_RANDMAX;
    return Z*gsl_ran_XXX_invRANDMAX;
}
inline double gsl_ran_XXX_uniform(void)
{
    return gsl_ran_XXX_uniform_wstate((void *)&state);
}

/* get/set randomState */
/* get/set randomState */
static inline void gsl_ran_XXX_copyState(void *tState,void *fState)
{
    bcopy((char *)fState,(char*)tState,sizeof(gsl_ran_XXX_randomState));
}

void *gsl_ran_XXX_getRandomState(void)
{
    gsl_ran_XXX_randomState *theState;
    theState = (gsl_ran_XXX_randomState *)calloc(1,sizeof(gsl_ran_XXX_randomState));
    gsl_ran_XXX_copyState(theState,(void *)&state);
    return (void *)theState;
}
void gsl_ran_XXX_setRandomState(void *vState)
{
    gsl_ran_XXX_copyState((void *)&state,(gsl_ran_XXX_randomState *)vState);
}


