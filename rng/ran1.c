/* The NR interface */
#include "gsl_ran.h"
float ran1(long *idum)
{
    if (*idum < 0) {
        *idum = -(*idum);
        gsl_ran_seed(*idum);
    }
    return gsl_ran_uniform();
}
