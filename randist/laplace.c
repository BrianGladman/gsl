/* $Id$ */
#include <math.h>
#include "gsl_ran.h"
#include "gsl_randist.h"

double gsl_ran_laplace(double mu)
{
    double u;
    do {
        u = 2 * gsl_ran_uniform() - 1.0;
    } while (u == 0.0);

    if (u < 0) {
      return  mu * log(-u) ;
    } else {
      return  -mu * log(u) ;
    }
}

/* The Laplace probability distribution is  
   
   p(x) = (1/(2 mu)) * exp( -|x/mu|)
   
   for -infty < x < infty  */

             

