/* $Id$ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_ran.h"
#include "gsl_randist.h"

int
main(int argc, char **argv)
{
    int i;
    double a=1;
    int n=1;
    double sum=0;
    int randseed=17;
    if (argc == 1) {
        printf("Usage: %s <a> <n> [seed]\n",argv[0]);
        printf("          Writes <n> random numbers\n");
        printf("          Gamma distributed with mean <a>,\n");
        printf("          optionally using <seed>\n");
        exit(0);
    }
    if (argc > 1) a = atof(argv[1]);
    if (argc > 2) n = atoi(argv[2]);
    if (argc > 3) {
        randseed = atoi(argv[3]);
        gsl_ran_seed(randseed);
    }

    for (i=0; i<n; ++i) {
        sum += gsl_ran_gamma(a);
    }
    printf("Average %g\n",sum/n);
    return 0;
}


