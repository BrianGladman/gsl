/* $Id$ */
/* Boilerplate: this is meant to be included in PNRG's */

void GSL_seed(int s)
{
    GSL_seed_wstate((void *)&state,s);
}
inline unsigned long GSL_random(void)
{
    return GSL_random_wstate((void *)&state);
}
inline double GSL_uniform(void)
{
    return GSL_uniform_wstate((void *)&state);
}

/* get/set randomState */
void *GSL_getRandomState(void)
{
    GSL_randomState *theState;
    theState = (GSL_randomState *)calloc(1,sizeof(GSL_randomState));
    GSL_copyRandomState(theState,(void *)&state);
    return (void *)theState;
}
void GSL_setRandomState(void *vState)
{
    GSL_copyRandomState((void *)&state,(GSL_randomState *)vState);
}


