/* $Id$ */
/* Boilerplate: this is meant to be included in PNRG's */

void gsl_ran_xxx_seed(int s)
{
    gsl_ran_xxx_seed_wstate((void *)&state,s);
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


