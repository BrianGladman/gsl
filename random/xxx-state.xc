/* $Id$ */
/* This file is generated automatically from the template file:
 * x.x.x-state.xc
 *
 * It defines functions gsl_ran_xxx_
 *   seed(),random(),uniform(),copyState(),getState(),setState()
 * in terms of the routines defined by the human programmer in xxx.c
 *
 * It is #include'd in file xxx.c
 */

void gsl_ran_xxx_seed(int s)
{
    gsl_ran_xxx_seed_wstate((void *)&state,s);
}
inline void gsl_ran_xxx_printState(void *s)
{
    gsl_ran_xxx_printState_p((gsl_ran_xxx_randomState *)s);
}
inline unsigned long gsl_ran_xxx_random(void)
{
    return gsl_ran_xxx_random_wstate((void *)&state);
}
inline double gsl_ran_xxx_uniform(void)
{
    return gsl_ran_xxx_uniform_wstate((void *)&state);
}

/* get/set randomState */
/* get/set randomState */
static inline void gsl_ran_xxx_copyState(void *tState,void *fState)
{
    bcopy((char *)fState,(char*)tState,sizeof(gsl_ran_xxx_randomState));
}

void *gsl_ran_xxx_getRandomState(void)
{
    gsl_ran_xxx_randomState *theState;
    theState = (gsl_ran_xxx_randomState *)calloc(1,sizeof(gsl_ran_xxx_randomState));
    gsl_ran_xxx_copyState(theState,(void *)&state);
    return (void *)theState;
}
void gsl_ran_xxx_setRandomState(void *vState)
{
    gsl_ran_xxx_copyState((void *)&state,(gsl_ran_xxx_randomState *)vState);
}


