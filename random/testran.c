/* $Id$ */
/* Test routine for the uniform random number generators */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_random.h"

int
main(int argc, char **argv)
{
    int i,n=1;
    unsigned long r,rmax;
    double sum,sigma;
    int randseed=17;
    void *tmpState;

    if (argc == 1) {
        printf("Usage: %s <n> [seed]\n",argv[0]);
        printf("          Tests random number generator\n");
        printf("          usiing <n> trials, and\n");
        printf("          optionally using <seed>\n");
        exit(0);
    }
    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) {
        /* these algorithms should work even if the seed is
         * not initially set! */
        randseed = atoi(argv[2]);
        GSL_seed(randseed);
    }

    printf("Nominal Maximum random  = %g\n",GSL_randmax());
    rmax=0;
    for (i=0; i<n; ++i) {
        r = GSL_random();
        if (rmax < r) rmax = r;
    }
    printf("Largest Observed random = %lu\n",rmax);
    printf("                        = %g\n",(double)rmax);
        
    sum=0;
    for (i=0; i<n; ++i)
        sum += GSL_uniform();
    sum /= n;
    /* expect sum to have variance == n*(1/12) */
    /* so average should have variance == 1/(12*n) */
    sigma = (sum - 0.5)*sqrt(12.0*n);
    printf("Sum test: %.2f sigmas\n",sigma);

    printf("Testing getRandomState/setRandomState:\n");
    printf("The following sets of numbers should be identical.\n");
    tmpState = GSL_getRandomState();
    for (i=0; i<5; ++i)
        printf("%c %.6f",(i==0 ? '(' : ','),GSL_uniform());
    printf(" )\n");
    GSL_setRandomState(tmpState); 
    for (i=0; i<5; ++i)
        printf("%c %.6f",(i==0 ? '(' : ','),GSL_uniform());
    printf(" )\n");
    GSL_setRandomState(tmpState);
    for (i=0; i<5; ++i) {
        printf("%c %.6f",(i==0 ? '(' : ','),GSL_uniform_wstate(tmpState));
        /* we have two parallel streams now, GSL_uniform() */
        GSL_uniform();
        GSL_uniform();
        GSL_uniform();
    }
    printf(" )\n");
    

    return 0;
}
    




