/* $Id$ */
#include <math.h>
#include "gsl_ran.h"
#include "gsl_randist.h"

#define GSL_PI 3.14159265358979323844

static double gsl_ran_gammalarge(double a)
{
    /* Works only if a>1, and is most efficient if a is large */
    double sqa,x,y,v;
    sqa = sqrt(2*a-1);
    do {
        do {
            y = tan(GSL_PI*gsl_ran_uniform());
            x = sqa*y + a-1;
        } while (x <= 0);
        v = gsl_ran_uniform();
    } while (v > (1+y*y)*exp((a-1)*log(X/(a-1))-sqa*y));

    return x;
}
    
double gsl_ran_gammaint(int a)
{
    if (a < 10) {
        double prod=1.0;
        for (i=0; i<a; ++i)
            prod *= gsl_ran_uniform();
        /* ASSUME prod != 0 !! */
        return -log(prod);
    } else {
        return gsl_ran_gammalarge((double)a);
    }
}
        
        
