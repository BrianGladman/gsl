/* $Id$ */
/* Benchmark routine for the uniform random number generators */
#include <stdio.h>
#include <stdlib.h>
#include "gsl_random.h"

int
main(int argc, char **argv)
{
    int i,n=1000000;
    int randseed=17;
    double sum;

    if (argc == 1) {
        printf("Usage: %s <n> [seed]\n",argv[0]);
        printf("          Computes sum of <n> random numbers\n");
        printf("          optionally using <seed>\n");
        exit(0);
    }
    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) randseed = atoi(argv[2]);

    GSL_seed(randseed);

    sum=0;
    for (i=0; i<n; ++i)
        sum += GSL_gauss();
    sum /= n;
    printf("Average of %d uniform random numbers: %.10f\n",n,sum);
    return 0;
}
    




