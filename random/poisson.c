/* $Id$ */
/* Get Poisson Deviates */

#include <math.h>
#include "gsl_ran.h"

double gsl_ran_exponential(double mu)
{
    double u;
    do {
        u = gsl_ran_uniform();
    } while (u<=0.0);
    return -mu*log(u);
}
    
    
    
int gsl_ran_poisson(double mu)
{
    /* This method works well when mu is small */
    double emu;
    double prod=1.0;
    int n=0;

    emu = exp(-mu);
    while (prod > emu) { 
        prod *= gsl_ran_uniform();
        ++n;
    }
    return n-1;
}

void gsl_ran_poisson_array(double mu, int N, int *p)
{
    int i,n;
    double emu;
    double prod=1.0;
    emu = exp(-mu);

    for (i=0; i<N; ++i) {
        n=0;
        while (prod > emu) { 
            prod *= gsl_ran_uniform();
            ++n;
        }
        p[i] = n-1;
    }
}    
