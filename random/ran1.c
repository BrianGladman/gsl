/* The NR interface */
#include "gsl_random.h"
float ran1(long *idum)
{
    if (*idum < 0) {
        *idum = -(*idum);
        GSL_seed(*idum);
    }
    return GSL_uniform();
}
