/* $Id$ */
/* Benchmark routine for Poisson random number generator */
#include <stdio.h>
#include <stdlib.h>
#include "gsl_ran.h"

int
main(int argc, char **argv)
{
    int i,n=1000000;
    int randseed=17;
    double mu=1.0;
    double sum=0.0;

    if (argc == 1) {
        printf("Usage: %s <mu> <n> [seed]\n",argv[0]);
        printf("          Computes the average of <n> random numbers\n");
        printf("          Poisson distributed wiht mean <mu>,\n");
        printf("          optionally using <seed>\n");
        exit(0);
    }
    if (argc > 1) mu = atof(argv[1]);
    if (argc > 2) n = atoi(argv[2]);
    if (argc > 3) randseed = atoi(argv[3]);

    gsl_ran_seed(randseed);

    sum=0;
    for (i=0; i<n; ++i)
        sum += gsl_ran_poisson(mu);
    sum /= n;
    printf("Average of %d Poisson random numbers: %.10f\n",n,sum);
    return 0;
}
    





