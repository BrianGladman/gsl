/* $Id$ */
/* Test routine for the gaussianian random number generators */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_ran.h"

int
main(int argc, char **argv)
{
    int i,n=1;
    unsigned long r,rmax;
    double sum,sigma;
    int randseed=17;
    gsl_ran_gaussianRandomState *tmpState;

    if (argc == 1) {
        printf("Usage: %s <n> [seed]\n",argv[0]);
        printf("          Tests random number generator\n");
        printf("          usiing <n> trials, and\n");
        printf("          optionally using <seed>\n");
        exit(0);
    }
    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) {
        randseed = atoi(argv[2]);
        gsl_ran_seed(randseed);
    }

    sum=0;
    for (i=0; i<n; ++i)
        sum += gsl_ran_gaussian();
    sum /= n;
    /* expect sum to have variance == n */
    /* so average should have variance == 1/n */
    sigma = sum*sqrt(n);
    printf("Sum test: %.2f sigmas\n",sigma);

    printf("Testing getRandomState/setRandomState:\n");
    printf("The following sets of numbers should be identical.\n");
    tmpState = gsl_ran_getGaussState();
    for (i=0; i<5; ++i)
        printf("%c %.6f",(i==0 ? '(' : ','),gsl_ran_gaussian());
    printf(" )\n");
    gsl_ran_setGaussState(tmpState); 
    for (i=0; i<5; ++i)
        /* gsl_ran_gaussian() doesn't influence tmpState */
        printf("%c %.6f",(i==0 ? '(' : ','),gsl_ran_gaussian());
    printf(" )\n");

    for (i=0; i<5; ++i) {
        printf("%c %.6f",(i==0 ? '(' : ','),gsl_ran_gaussian_wstate(tmpState));
        /* these should not influence tmpState */
        gsl_ran_uniform();
        gsl_ran_gaussian();
    }
    printf(" )\n");
    

    return 0;
}
    




