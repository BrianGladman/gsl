/* $Id$ */
#include <math.h>
#include "gsl_ran.h"
#include "gsl_randist.h"

double gsl_ran_exponential(double mu)
{
    double u;
    do {
        u = gsl_ran_uniform();
    } while (u<=0.0);
    return -mu*log(u);
}
    
    
